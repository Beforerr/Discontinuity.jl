process!(df::AbstractDataFrame) = begin
    df |>
    dropmissing |> # TODO: keep missing values
    keep_good_fit! |>
    standardize_df! |>
    compute_params! |>
    compute_Alfvenicity_params!
end

"""the MVAB method can achieve acceptable accuracy when either |B|/|B| > 0.05 or ω > 60°. @liuFailuresMinimumVariance2023"""
assign_mva_quality!(df) = @rtransform!(df, :mva_quality = (:ω > 60) | (:dBmag_over_Bmag > 0.05))
function filter_low_mva_quality!(df)
    "mva_quality" in names(df) || assign_mva_quality!(df)
    filter!(:mva_quality => ==(true), df)
end

function compute_orientation_params!(df)
    @chain df begin
        assign_mva_quality!
        @rtransform! begin
            :θ_mva_cross = angle_between_90(:n_mva, :n_cross)
            :B_n_mva_norm = abs(:B_n_mva / ustrip(:B_mag))
            :B_n_cross_norm = abs(:B_n_cross / ustrip(:B_mag))
        end
    end
end

function compute_params!(df; l_unit=DEFAULT_L_UNIT, j_unit=DEFAULT_J_UNIT)
    cols = names(df)
    "dB" in cols && @chain df begin
        @rtransform!(
            :ω = vector_angle(:"B.vec.before", :"B.vec.after"),
            :ω_in = vector_angle(:"B.vec.before"[1:2], :"B.vec.after"[1:2]),
        )
        @transform!(
            :"B.mean" = (:"B.before" .+ :"B.after") ./ 2,
            :"n.mean" = (:"n.before" .+ :"n.after") ./ 2,
        )
        @transform!(
            :dB_over_B = ((:"B.before" .- :"B.after") ./ :"B.mean"),
            :dn_over_n = ((:"n.before" .- :"n.after") ./ :"n.mean"),
        )
    end

    "V" in cols && @transform! df @astable begin
        :V_mag = norm.(:V)
        :V_n_cross = sproj.(:V, :n_cross)
        :L_n_cross = @. abs(:duration * :V_n_cross) |> l_unit
        :J_m_max_cross = @. abs(gradient_current(:grad, :V_n_cross)) |> j_unit
        :V_n_mva = sproj.(:V, :n_mva)
        :L_n_mva = @. abs(:duration * :V_n_mva) |> l_unit
        :J_m_max_mva = @. abs(gradient_current(:grad, :V_n_mva)) |> j_unit
    end

    "n" in cols && @transform! df @astable begin
        :n = _unitify_n.(:n)
        :d_i = @. inertial_length(:n) |> l_unit
        :V_A = Alfven_speed.(:B_mag, :n) # Alfven speed
        :J_A = @. upreferred(:V_A * :n * Unitful.q)

        :L_n_cross_norm = @. abs(:L_n_cross / :d_i)
        :L_n_mva_norm = @. abs(:L_n_mva / :d_i)

        :V_A_lmn_before = alfven_velocity.(:B_lmn_before, :n)
        :V_A_lmn_after = alfven_velocity.(:B_lmn_after, :n)
        :J_m_max_mva_norm = NoUnits.(:J_m_max_mva ./ :J_A)
        :J_m_max_cross_norm = NoUnits.(:J_m_max_cross ./ :J_A)
    end

    "T" in cols && @transform! df @astable begin
        :β = plasma_beta.(:T, :n, :B_mag)
    end
    "T.before" in cols && @transform! df :"T.mean" = (:"T.before" .+ :"T.after") ./ 2
    return df
end

function compute_Alfvenicity_params!(df)
    @chain df begin
        @rtransform!(
            :V_us_l = sproj(:V_us, :e_max),
            :V_ds_l = sproj(:V_ds, :e_max),
            :B_us_l = sproj(:"B.vec.before", :e_max),
            :B_ds_l = sproj(:"B.vec.after", :e_max)
        )
        @transform!(
            :V_A_us_l = Alfven_velocity.(:B_us_l, :"n.before"),
            :V_A_ds_l = Alfven_velocity.(:B_ds_l, :"n.after"),
        )
        @transform!(
            :"v.Alfven.change.l" = abs.(:V_A_us_l .- :V_A_ds_l),
            :"v.ion.change.l" = abs.(:V_us_l .- :V_ds_l)
        )
        @transform! :v_l_ratio = :"v.ion.change.l" ./ :"v.Alfven.change.l"
        @transform! :Λ_t = 1 .- :v_l_ratio .^ 2
    end

    if "v.Alfven.change.l.fit" in names(df)
        @transform! df :v_l_fit_ratio = :"v.ion.change.l" ./ :"v.Alfven.change.l.fit"
    end
    return df
end