"""Calculate Sonnerup (2018) quality index for Walén relation assessment."""
function sonnerup_quality_index(ΔV, ΔV_A)
    denominator = norm(ΔV) + norm(ΔV_A)
    # Q⁺: ΔV - ΔV_A (co-aligned case)
    Q⁺ = 1 - norm(ΔV - ΔV_A) / denominator
    # Q⁻: ΔV + ΔV_A (anti-aligned case)
    Q⁻ = -(1 - norm(ΔV + ΔV_A) / denominator)
    return abs(Q⁺) > abs(Q⁻) ? Q⁺ : Q⁻
end

function compute_Alfvenicity_maximum_variance_direction!(df)
    return @transform! df @astable begin
        ΔV_l = passmissing(sproj).(:ΔV, :e_max)
        ΔV_A_l = passmissing(sproj).(:ΔV_A, :e_max)
        :V_l_ratio = @. abs(ΔV_l / ΔV_A_l) |> NoUnits
        :V_l_ratio_max = @. abs(:ΔV_l_max / ΔV_A_l) |> NoUnits
    end
end

function compute_Alfvenicity_params!(df)
    @chain df begin
        @transform! @astable begin
            # Calculate full velocity and Alfvén velocity changes
            :ΔV = passmissing(-).(:V_ds, :V_us)
            :ΔV_A = Alfven_velocity.(:B_ds, :n) .- Alfven_velocity.(:B_us, :n)
            :ΔV_ratio = @. norm(:ΔV) / norm(:ΔV_A) |> NoUnits
            :ΔV_cosθ = passmissing(cosd ∘ angle_between).(:ΔV, :ΔV_A)
            :Q_sonnerup = passmissing(sonnerup_quality_index).(:ΔV, :ΔV_A)
        end
        compute_Alfvenicity_maximum_variance_direction!(_)
        @transform! @astable begin
            :Λ_t = 1 .- :V_l_ratio .^ 2
        end
    end

    if "v.Alfven.change.l.fit" in names(df)
        @transform! df :v_l_fit_ratio = :"v.ion.change.l" ./ :"v.Alfven.change.l.fit"
    end
    return df
end
