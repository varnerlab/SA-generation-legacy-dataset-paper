# ──────────────────────────────────────────────────────────────────────────────
# test_revision2_helpers.jl
#
# Tests for revision2_helpers.jl (npj round-2 minor revision, Reviewer 3).
# The first two testsets pin the helpers to numbers already published:
#   - the exact-sampler reference row (exact_sampler_reference.csv), and
#   - the 54/100 biological pass count of the published cohort.
# The rest check the new constraint, rejection, Langevin, and bandwidth code.
#
# Run from code/:  julia experiments/test_revision2_helpers.jl
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests, Test
include(joinpath(@__DIR__, "mechanistic_eval.jl"))
include(joinpath(@__DIR__, "revision2_helpers.jl"))

const ref = load_reference()   # once, at top level; every testset below uses it

@testset "revision2 helpers" begin

@testset "metric copy reproduces exact_sampler_reference.csv" begin
    dirs, _ = exact_directions(ref, ref.β_star, MersenneTwister(42), 100)
    df = decode_cohort(ref, dirs, ref.shared_magnitudes)
    m = evaluate_cohort(ref, dirs, df; with_mech=true)
    row = CSV.read(joinpath(_PATH_TO_DATA, "exact_sampler_reference.csv"), DataFrame)[1, :]
    for f in (:Median_MRE, :CrossVisit_Frob, :Mech_Overlap_TFonly, :Median_Novelty,
              :DCR_synth_to_real_median, :DCR_real_to_real_median)
        @test isapprox(getproperty(m, f), row[f]; rtol=1e-10)
    end
end

@testset "bio check reproduces 54/100 on the published cohort" begin
    @test count_bio_pass_published(ref) == 54
end

@testset "hormone affine map matches decode_sample" begin
    A = hormone_affine(ref)
    @test length(A.idx) == 6
    @test size(A.c) == (6, size(ref.X̂, 1))
    @test A.s == ref.std_params_concat.σ[A.idx]
    rng = MersenneTwister(7)
    for _ in 1:20
        u = normalize(randn(rng, size(ref.X̂, 1)))
        r = rand(rng, ref.pca_norms)
        full = decode_sample(r .* u, ref.pca_concat, ref.std_params_concat)
        @test isapprox(full[A.idx], A.a .+ r .* (A.c * u); rtol=1e-10, atol=1e-10)
    end
end

@testset "penalty is zero when hormones are nonnegative; gradient matches finite differences" begin
    A = hormone_affine(ref)
    rng = MersenneTwister(11)
    rmax = maximum(ref.pca_norms)
    n_zero = 0
    n_pos = 0
    for _ in 1:400
        ξ = ref.X̂[:, rand(rng, 1:size(ref.X̂, 2))] .+ 0.6 .* randn(rng, size(ref.X̂, 1))
        x = A.a .+ rmax .* (A.c * normalize(ξ))
        Cval = positivity_penalty(ξ, rmax, A)
        if all(x .>= 0)
            @test Cval == 0.0
            n_zero += 1
        else
            @test Cval > 0
            n_pos += 1
            g = positivity_grad(ξ, rmax, A)
            h = 1e-6
            g_fd = similar(ξ)
            for i in eachindex(ξ)
                e = zeros(length(ξ)); e[i] = h
                g_fd[i] = (positivity_penalty(ξ .+ e, rmax, A) - positivity_penalty(ξ .- e, rmax, A)) / (2h)
            end
            @test isapprox(g, g_fd; rtol=1e-5, atol=1e-8)
        end
    end
    @test n_zero > 0      # both branches were exercised
    @test n_pos > 0
end

@testset "hard rejection (λ = Inf) gives nonnegative hormones" begin
    A = hormone_affine(ref)
    rng = MersenneTwister(13)
    for _ in 1:200
        r = rand(rng, ref.pca_norms)
        ξ, tries = rejection_direction(ref.X̂, r, A; β=ref.β_star, λ=Inf, rng=rng)
        @test tries >= 1
        full = decode_sample(r .* normalize(ξ), ref.pca_concat, ref.std_params_concat)
        @test all(full[A.idx] .>= 0)
    end
    @test_throws ErrorException rejection_direction(ref.X̂, 1.0, A; β=ref.β_star, λ=Inf,
                                                    rng=MersenneTwister(1), max_tries=0)
end

@testset "leave-one-out bandwidth is an interior maximum; Scott formula" begin
    bw = kde_bandwidths(ref.X̂)
    d, K = size(ref.X̂)
    @test 0.01 < bw.loo < 3.0
    @test loo_loglik(ref.X̂, bw.loo) >= loo_loglik(ref.X̂, 0.9 * bw.loo)
    @test loo_loglik(ref.X̂, bw.loo) >= loo_loglik(ref.X̂, 1.1 * bw.loo)
    @test isapprox(bw.scott, bw.sigma * K^(-1 / (d + 4)); rtol=1e-12)
    @test isapprox(bw.silverman, bw.sigma * (4 / (d + 2))^(1 / (d + 4)) * K^(-1 / (d + 4)); rtol=1e-12)
end

@testset "untransformed representation reproduces the canonical memory" begin
    rep = build_representation(ref, Int[])
    @test size(rep.X̂) == size(ref.X̂)
    @test isapprox(rep.X̂, ref.X̂; rtol=1e-10, atol=1e-10)
    @test isapprox(rep.pca_norms, ref.pca_norms; rtol=1e-10)
    @test isapprox(rep.β_star, ref.β_star; rtol=1e-10)
    rng = MersenneTwister(23)
    for _ in 1:10
        y = randn(rng, size(ref.X̂, 1)) .* 10
        @test isapprox(decode_rep(rep, y), decode_sample(y, ref.pca_concat, ref.std_params_concat); rtol=1e-10, atol=1e-8)
    end
end

@testset "plausibility summary reproduces the published counts" begin
    s = plausibility_summary(ref, published_cohort_long(ref))
    @test s.Bio_Pass == 54
    @test s.Hormone_Neg_Profiles == 46
    @test s.Strict_Mech_Pass == 41
    @test s.Both_Pass == 24
end

@testset "log-hormone representation decodes to positive hormones" begin
    A = hormone_affine(ref)
    rep = build_representation(ref, A.idx)
    rng = MersenneTwister(29)
    for _ in 1:200
        u = normalize(randn(rng, size(rep.X̂, 1)))
        x = decode_rep(rep, rand(rng, rep.pca_norms) .* u)
        @test all(x[A.idx] .> 0)
    end
end

@testset "λ = 0 constrained ULA matches the exact mixture" begin
    A = hormone_affine(ref)
    n = 1000
    rng = MersenneTwister(17)
    d, K = size(ref.X̂)
    ula = Matrix{Float64}(undef, d, n)
    for i in 1:n
        k0 = rand(rng, 1:K)
        ξ0 = normalize(ref.X̂[:, k0] .+ 0.01 .* randn(rng, d))
        ula[:, i] = normalize(constrained_ula(ref.X̂, ξ0, 1.0, A; β=ref.β_star, λ=0.0,
                                              α=0.01, T=2000, rng=rng))
    end
    ex = exact_sample(ref.X̂, n; β=ref.β_star, rng=MersenneTwister(19)).Ξ
    exu = reduce(hcat, [normalize(ex[i, :]) for i in 1:n])
    for j in 1:2
        @test pvalue(ApproximateTwoSampleKSTest(ula[j, :], exu[j, :])) > 0.01
    end
end

end
