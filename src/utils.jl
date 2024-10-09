using LinearAlgebra

"""
    keep_good_fit!(df; rsquared=0.9)

Keep only the rows with a fit statistic R² greater than the given value.
"""
function keep_good_fit!(df; rsquared=0.9)
    @subset!(df, :"fit.stat.rsquared" .> rsquared)
end

function compute_params!(df)
    df = @chain df begin
        @rtransform!(
            :dB_norm = norm([:dB_x, :dB_y, :dB_z])
        )
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
        )
    end

    if "T.before" in names(df)
        @transform! df :"T.mean" = (:"T.before" .+ :"T.after") ./ 2
    end
    return df
end

function compute_Alfvenicity_params!(df)
    @chain df begin
        @transform!(
            :"v.Alfven.change.l" = abs.(:"v.Alfven.change.l"),
            :"v.ion.change.l" = abs.(:"v.ion.change.l")
        )
        @transform! :v_l_ratio = :"v.ion.change.l" ./ :"v.Alfven.change.l"
        @transform! :Λ_t = 1 .- :v_l_ratio .^ 2
    end

    if "v.Alfven.change.l.fit" in names(df)
        @transform! df :v_l_fit_ratio = :"v.ion.change.l" ./ :"v.Alfven.change.l.fit"
    end
    return df
end