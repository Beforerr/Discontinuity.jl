export var_mapping

# %%
# Define the labels for the plots
j_lab = L"Current Density ($nA/m^2$)"
l_lab = L"Thickness ($km$)"

l_norm_lab = "Normalized Thickness"
j_norm_lab = "Normalized Current Density"

di_lab = L"Ion Inertial Length ($km$)"
jA_lab = L"Alfvénic Current Density ($nA/m^2$)"

B_lab = L"Magnetic Field ($nT$)"
b_fit_lab = L"Fitted Amplitude ($nT$)"
density_lab = L"Density ($cm^{-3}$)"

# %%
# Define the mappings
di_map = :ion_inertial_length => di_lab
di_log_map = :ion_inertial_length => log10 => L"Log %$di_lab"

jA_map = :j_Alfven => jA_lab
jA_log_map = :j_Alfven => log10 => L"Log %$jA_lab"

v_Alfven_map = "v.Alfven.change.l" => L"\Delta V_{A,l}"
v_ion_map = "v.ion.change.l" => L"\Delta V_{i,l}"
v_l_ratio_map = "v_l_ratio" => L"\Delta V_{i,l} / \Delta V_{A,l}"

function var_mapping(; log_str="Log ")

    l_map = (;
        l=:L_n_cross => l_lab,
        l_log=:L_n_cross => log10 => LaTeXString("$log_str$l_lab"),
        l_norm=:L_n_cross_norm => l_norm_lab,
        l_norm_log=:L_n_cross_norm => log10 => LaTeXString("$log_str$l_norm_lab"),
    )

    j_map = (;
        j=:j0_k => j_lab,
        j_log=:j0_k => log10 => LaTeXString("$log_str$j_lab"),
        j_norm=:j0_k_norm => j_norm_lab,
        j_norm_log=:j0_k_norm => log10 => LaTeXString("$log_str$j_norm_lab"),
    )

    param_map = (;
        dB_over_B=:dB_over_B => L"\Delta B / B",
        dB_norm_over_B=("dB_norm", "B.mean") => (/) => L"|\Delta \mathbf{B} |/B",
        ω=:ω => "rotation angle",
        ω_in=:ω_in => "in-plane rotation angle",
        bn=:bn_over_b => abs => L"B_N/B",
        jA=:j_Alfven => jA_lab,
        jA_log=:j_Alfven => log10 => LaTeXString("$log_str$jA_lab"),

        # Parameters
        density=:"n.mean" => density_lab,
        density_log=:"n.mean" => log10 => LaTeXString("$log_str$density_lab"),
        B=:"B.mean" => B_lab,
        B_log=:"B.mean" => log10 => LaTeXString("$log_str$B_lab"),
        bm0_log=:"fit.vars.amplitude" => log10 ∘ abs => LaTeXString("$log_str$b_fit_lab"),
        beta=:β => log10 => L"%$(log_str)Plasma Beta $\beta$"
    )
    return merge(l_map, j_map, param_map)
end

# baremodule DefaultMapping
# all_maps = ...
# for m in all_maps
#     Core.eval(DefaultMapping, Expr(:import, Expr(:(.), :Discontinuity, m)))
#     Core.eval(DefaultMapping, Expr(:export, m))
# end
# end