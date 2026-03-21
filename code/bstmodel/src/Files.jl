function write_legacy_coagulation_population_file(patients::Dict{Int64,Array{Float64,1}}, path::String; visitid = 1)

    # initialize -
    df = DataFrame(
        ssubjid = Int64[],
        visitid = Int64[],
        II = Float64[],
        V = Float64[],
        VII = Float64[],
        VIII = Float64[],
        IX = Float64[],
        X = Float64[],
        XI = Float64[],
        XII = Float64[],
        AT = Float64[],
        PC = Float64[],
        TFPI = Float64[],
        Lagtime = Float64[],
        Peak = Float64[],
        TPeak = Float64[],
        Max = Float64[],
        EPT = Float64[],
    )

    # process each synthetic patient -
    for (i,data) ∈ patients

        # grab -
        results_tuple = (
            ssubjid = i,
            visitid = visitid,
            II = data[1],
            V = data[2],
            VII = data[3],
            VIII = data[4],
            IX = data[5],
            X = data[6],
            XI = data[7],
            XII = data[8],
            AT = data[9],
            PC = data[10],
            TFPI = data[11],
            Lagtime = data[12],
            Peak = data[13],
            TPeak = data[14],
            Max = data[15],
            EPT = data[16]
        );
        
        push!(df, results_tuple)
    end

    # save -
    CSV.write(path, df)
end