# %%
# Define the mappings
di_map = :ion_inertial_length => di_lab
di_log_map = :ion_inertial_length => log10 => L"Log %$di_lab"
l_map = :L_k => l_lab
l_norm_map = :L_k_norm => l_norm_lab
l_log_map = :L_k => log10 => L"Log %$l_lab"
l_norm_log_map = :L_k_norm => log10 => L"Log %$l_norm_lab"

jA_map = :j_Alfven => jA_lab
jA_log_map = :j_Alfven => log10 => L"Log %$jA_lab"

j_map = :j0_k => j_lab
j_log_map = :j0_k => log10 => L"Log %$j_lab"
j_norm_map = :j0_k_norm => j_norm_lab
j_norm_log_map = :j0_k_norm => log10 => L"Log %$j_norm_lab"

v_Alfven_map = "v.Alfven.change.l" => L"\Delta V_{A,l}"
v_ion_map = "v.ion.change.l" => L"\Delta V_{i,l}"