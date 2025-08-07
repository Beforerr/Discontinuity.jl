include("processing/Alfvenicity.jl")

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
function filter_low_mva_quality(df; view = true)
    "mva_quality" in names(df) || assign_mva_quality!(df)
    return filter(:mva_quality => identity, df; view)
end

function compute_orientation_params!(df)
    return @chain df begin
        assign_mva_quality!
        @rtransform! begin
            :θ_mva_cross = angle_between_90(:n_mva, :n_cross)
            :B_n_mva_norm = abs(:B_n_mva / ustrip(:B_mag))
            :B_n_cross_norm = abs(:B_n_cross / ustrip(:B_mag))
        end
    end
end

function compute_params!(df; l_unit = L_UNIT, j_unit = J_UNIT)
    cols = names(df)

    @rtransform! df @astable begin
        :ω_in = angle_between(:B_lmn_after[1:2], :B_lmn_before[1:2])
    end

    "V" in cols && @transform! df @astable begin
        :V_mag = norm.(:V)
        :V_n_cross = sproj.(:V, :n_cross)
        :L_n_cross = @. abs(:duration * :V_n_cross) |> l_unit
        :J_m_max_cross = @. abs(gradient_current(:grad, :V_n_cross)) |> j_unit
        :V_n_mva = sproj.(:V, :n_mva)
        :V_A_n_mva = Alfven_velocity.(:B_n_mva, :n)
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

        :J_m_max_mva_norm = NoUnits.(:J_m_max_mva ./ :J_A)
        :J_m_max_cross_norm = NoUnits.(:J_m_max_cross ./ :J_A)
    end

    "T" in cols && @transform! df @astable begin
        :β = plasma_beta.(:T, :n, :B_mag)
    end
    "T.before" in cols && @transform! df :"T.mean" = (:"T.before" .+ :"T.after") ./ 2
    return assign_mva_quality!(df)
end

skipmissing_sum(args...) = all(ismissing, args) ? missing : sum(skipmissing(args))

calc_v_l_ratio_Λ(v_l_ratio::T, Λ::Number) where {T} = (Λ > 1) ? T(NaN) : v_l_ratio / sqrt(1 - Λ)
calc_v_l_ratio_Λ(v_l_ratio, Λ::Missing) = missing

"""
compute_anisotropy_params!(df, :ion => (:Tp_para, :Tp_perp), :electron => (:Te_para, :Te_perp))
"""
function compute_anisotropy_params!(df, T3s::Pair...)
    Λ_ss = map(T3s) do (k, v)
        Tpara = df[!, v[1]]
        Tperp = df[!, v[2]]
        Λ_sym = Symbol(:Λ_, k)
        Λ_s = anisotropy.(df.B_mag, df.n, Tpara, Tperp)
        df[!, Λ_sym] = Λ_s
    end
    Λ_all = map(skipmissing_sum, Λ_ss...)
    if "V_l_ratio" in names(df)
        @transform! df :V_l_ratio_Λ = calc_v_l_ratio_Λ.(:V_l_ratio, Λ_all)
    end
    df.Λ = Λ_all
    return df
end

function compute_params_py!(df)
    return "dB" in cols && @chain df begin
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
        @rtransform!(
            :V_us_l = sproj(:V_us, :e_max),
            :V_ds_l = sproj(:V_ds, :e_max),
            :B_us_l = sproj(:"B.vec.before", :e_max),
            :B_ds_l = sproj(:"B.vec.after", :e_max)
        )
    end
end


function waiting_time(time; δt=Dates.Minute(1))
    # unique and order the time
    τ = diff(time |> unique |> sort)
    return τ ./ δt
end

waiting_time(df::AbstractDataFrame, col=:time; kwargs...) = waiting_time(df[!, col]; kwargs...)