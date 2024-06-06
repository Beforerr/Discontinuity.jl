
using FileIO

function keep_good_fit!(df; rsquared=0.9)
    @subset!(df, :"fit.stat.rsquared" .> rsquared)
end

function process!(df::AbstractDataFrame)
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
        @transform! :v_l_ratio = :"v.ion.change.l" ./ :"v.Alfven.change.l"
        @transform! :Λ_t = 1 .- :v_l_ratio .^ 2
        unique!(["t.d_start", "t.d_end"])
    end

    if "T.before" in names(df)
        @transform! df :"T.mean" = (:"T.before" .+ :"T.after") ./ 2
    end

    df |> keep_good_fit!
end

function load(path; )
    df = path |> Arrow.Table |> DataFrame |> dropmissing
    df |> process!
end