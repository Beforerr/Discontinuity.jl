using Arrow
using DataFrames,
    DataFramesMeta

function keep_good_fit(df)
    df = @chain df begin
        filter(:"fit.stat.rsquared" => >(0.9), _)
    end
end

function process(df)
    df = @chain df begin
        subset!(names(df) .=> ByRow(isfinite))
        transform!(names(df, Float32) .=> ByRow(Float64); renamecols=false) # Convert all columns of Float32 to Float64
        @transform!(
            :"B.mean" = (:"B.before" .+ :"B.after") ./ 2,
            :"n.mean" = (:"n.before" .+ :"n.after") ./ 2,
        )
        @transform!(
            :dB_over_B = (:"B.change" ./ :"B.mean"),
            :dn_over_n = (:"n.change" ./ :"n.mean"),
        )
        @transform!(
            :j0_k = abs.(:j0_k),
            :j0_k_norm = abs.(:j0_k_norm),
            :"v.Alfven.change.l" = abs.(:"v.Alfven.change.l"),
            :"v.ion.change.l" = abs.(:"v.ion.change.l")
        )
        @transform! :Λ_t = 1 .- (:"v.ion.change.l" ./ :"v.Alfven.change.l") .^ 2
        unique!(["t.d_start", "t.d_end"])
    end

    if "T.before" in names(df)
        @transform! df :"T.mean" = (:"T.before" .+ :"T.after") ./ 2
    end

    df |> keep_good_fit
end

function load(path::String)
    df = path |> Arrow.Table |> DataFrame |> dropmissing
    df |> process
end