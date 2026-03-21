function _loss_scalar(κ::Array{Float64,1}, Y::Array{Float64,1},  model::BSTModel, mvar::Union{Nothing,Array{Float64,1}})

    # get the G matrix -
    G = model.G

    # κ map
    # 1 - 9 : α vector
    model.α = κ[1:9]
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = κ[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = κ[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = -1*κ[12]
    model.G = G

    # run the model -
    (T,Xₘ) = evaluate(model, tspan=(0.0, 40.0));

    # compute the model ouput -
    Yₘ = model_output_vector(T,Xₘ[:,9]) # properties of the FIIa curve

    # compute the fitness -
    ϵ = nothing;
    if (isnothing(mvar) == false)

        # compute the scaled squared error -
        number_of_outputs = length(Yₘ);
        tmp = zeros(number_of_outputs);
        for i ∈ 1:number_of_outputs
            tmp[i] = (Y[i] - Yₘ[i])^2 / mvar[i]
        end
        ϵ = sum(tmp);
    else
        # normalized relative squared error — each feature contributes equally
        number_of_outputs = length(Yₘ);
        tmp = zeros(number_of_outputs);
        for i ∈ 1:number_of_outputs
            if abs(Y[i]) > 1e-12
                tmp[i] = ((Y[i] - Yₘ[i]) / Y[i])^2
            else
                tmp[i] = (Y[i] - Yₘ[i])^2
            end
        end
        ϵ = sum(tmp);
    end

    # return -
    return ϵ
end

# ──────────────────────────────────────────────────────────────────────────────
# _setup_model! — configure static factors and initial conditions from data
# ──────────────────────────────────────────────────────────────────────────────
function _setup_model!(model::BSTModel, index::Int, training_df::DataFrame)

    SF = 1e9

    # setup static factors
    sfa = model.static_factors_array
    sfa[1] = training_df[index,:TFPI]       # 1 TFPI
    sfa[2] = training_df[index,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF               # 3 TF
    sfa[6] = 0.005                      # 6 TRAUMA

    # setup initial conditions from patient data
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[index, :II]         # 1 FII
    xₒ[2] = training_df[index, :VII]        # 2 FVII
    xₒ[3] = training_df[index, :V]          # 3 FV
    xₒ[4] = training_df[index, :X]          # 4 FX
    xₒ[5] = training_df[index, :VIII]       # 5 FVIII
    xₒ[6] = training_df[index, :IX]         # 6 FIX
    xₒ[7] = training_df[index, :XI]         # 7 FXI
    xₒ[8] = training_df[index, :XII]        # 8 FXII
    xₒ[9] = (1e-14)*SF                      # 9 FIIa
    model.initial_condition_array = xₒ

    return nothing
end

# ──────────────────────────────────────────────────────────────────────────────
# _apply_parameters! — map parameter vector κ into the BST model
# ──────────────────────────────────────────────────────────────────────────────
function _apply_parameters!(model::BSTModel, κ::Array{Float64,1})

    model.α = κ[1:9]
    G = model.G
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = κ[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = κ[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = -1*κ[12]
    model.G = G

    return nothing
end

# ──────────────────────────────────────────────────────────────────────────────
# Default parameter bounds (columns: default, lower, upper)
# ──────────────────────────────────────────────────────────────────────────────
const _DEFAULT_PARAMETER_BOUNDS = [
    # default   lower   upper
    0.061   0.0001 25.0   ; # 1  α₁ (FVII → FVIIa)
    1.0     0.0001 25.0   ; # 2  α₂ (FVIII → FVIIIa)
    1.0     0.0001 25.0   ; # 3  α₃ (FIX → FIXa)
    1.0     0.0001 25.0   ; # 4  α₄ (FV → FVa)
    1.0     0.0001 25.0   ; # 5  α₅ (FX → FXa)
    1.0     0.0001 25.0   ; # 6  α₆ (FII → FIIa)
    1.0     0.0001 25.0   ; # 7  α₇ (FXI → FXIa)
    1.0     0.0001 25.0   ; # 8  α₈ (FXII → FXIIa)
    0.70    0.0001 25.0   ; # 9  α₉ (FIIa → FIIa_inactive)
    0.11    0.0001 25.0   ; # 10 G: TF effect on FVIIa
    0.045   0.0001 25.0   ; # 11 G: AT deactivation
    0.065   0.0001 25.0   ; # 12 G: TFPI inhibition
]

# ──────────────────────────────────────────────────────────────────────────────
# learn_bst_model_parameters — multi-start optimization with normalized loss
#
# Runs n_restarts independent Nelder-Mead optimizations from random starting
# points (log-uniform within bounds), keeps the best result.
# ──────────────────────────────────────────────────────────────────────────────
function learn_bst_model_parameters(index::Int, model::BSTModel, training_df::DataFrame;
    pₒ::Union{Nothing,Array{Float64,1}} = nothing,
    mvar::Union{Nothing,Array{Float64,1}} = nothing,
    n_restarts::Int = 10,
    iterations_per_restart::Int = 1000,
    time_limit_per_restart::Float64 = 300.0,
    show_trace::Bool = false)

    # setup model with patient-specific ICs and static factors
    _setup_model!(model, index, training_df)

    # target output vector
    Y = Array{Float64,1}(undef,5)
    Y[1] = training_df[index, :Lagtime]
    Y[2] = training_df[index, :Peak]
    Y[3] = training_df[index, :TPeak]
    Y[4] = training_df[index, :Max]
    Y[5] = training_df[index, :EPT]

    # parameter bounds
    κ = copy(_DEFAULT_PARAMETER_BOUNDS)
    P = size(κ, 1)
    lb = κ[:, 2]
    ub = κ[:, 3]

    # objective function (uses normalized relative error by default)
    obj_function(p) = _loss_scalar(p, Y, model, mvar)

    # track best across restarts
    best_p = nothing
    best_ϵ = Inf

    inner_optimizer = NelderMead()

    for restart in 1:n_restarts

        # generate starting point
        if restart == 1 && !isnothing(pₒ)
            # first restart: use provided starting point
            p_start = copy(pₒ)
            for i in 1:P
                p_start[i] = clamp(p_start[i], lb[i], ub[i])
            end
        elseif restart == 1
            # first restart: use defaults with small perturbation
            p_start = κ[:, 1] .* (1.0 .+ 0.1 .* randn(P))
            for i in 1:P
                p_start[i] = clamp(p_start[i], lb[i], ub[i])
            end
        else
            # subsequent restarts: log-uniform random within bounds
            p_start = zeros(P)
            for i in 1:P
                log_lb = log10(lb[i])
                log_ub = log10(ub[i])
                p_start[i] = 10.0^(log_lb + rand() * (log_ub - log_lb))
            end
        end

        try
            results = optimize(obj_function, lb, ub, p_start, Fminbox(inner_optimizer),
                Optim.Options(time_limit=time_limit_per_restart, show_trace=show_trace,
                              show_every=100, iterations=iterations_per_restart))

            ϵ = Optim.minimum(results)
            if ϵ < best_ϵ
                best_ϵ = ϵ
                best_p = Optim.minimizer(results)
            end
        catch e
            # silently skip failed restarts (ODE solver divergence, etc.)
            continue
        end
    end

    # if no restart succeeded, error
    if best_p === nothing
        error("All $n_restarts optimization restarts failed for patient $index")
    end

    # run final simulation with best parameters
    _apply_parameters!(model, best_p)
    _setup_model!(model, index, training_df)  # re-setup ICs (optimizer may have mutated)
    (T, Xₘ) = evaluate(model, tspan=(0.0, 40.0))
    Yₘ = model_output_vector(T, Xₘ[:, 9])

    return (best_p, T, Xₘ, Yₘ, Y, best_ϵ)
end
