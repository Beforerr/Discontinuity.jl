process!(df::AbstractDataFrame) = begin
    df |>
    dropmissing |> # TODO: keep missing values
    keep_good_fit! |>
    standardize_df! |>
    compute_params! |>
    compute_Alfvenicity_params!
end

function compute_params!(df)
    df = @chain df begin
        @rtransform!(
            :dB_norm = norm(:dB),
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
        @transform!(
            :j0_k = abs.(:j0_k),
            :j0_k_norm = abs.(:j0_k ./ :j_Alfven),
        )
    end

    :"n_cross" in names(df) && @rtransform!(
        df,
        :θ_nk = vector_angle(:e_min, :n_cross),
        :L_n_cross_norm = abs(:L_n_cross / :ion_inertial_length),
    )

    if "T.before" in names(df)
        @transform! df :"T.mean" = (:"T.before" .+ :"T.after") ./ 2
    end
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