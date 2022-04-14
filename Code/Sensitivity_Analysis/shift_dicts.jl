include("../../Data/exp_data.jl")

tvals_dict = Dict("Snf1" => t_Snf1, "cAMP" => t_cAMP, "sch9Delta_cAMP" => t_sch9Delta_cAMP,
"Sch9P_glutamine_H" => t_Sch9P_glutamine_H, "Sch9P_glutamine_L" => t_Sch9P_glutamine_L,
"Sch9_gtr1Delta" => t_Sch9_gtr1Delta, "Sch9_glucose_starve" => t_Sch9_glucose_starve,
"Sch9_glucose_relief" => t_Sch9_glucose_relief, "Mig1_glucose_relief" => t_Mig1_glucose_relief,
"Rib_rap" => t_Rib_rap)

dat_dict = Dict("Snf1" => raw_data_Snf1, "cAMP" => raw_data_cAMP, "sch9Delta_cAMP" => raw_data_sch9Delta_cAMP,
"Sch9P_glutamine_H" => raw_data_Sch9P_glutamine_H, "Sch9P_glutamine_L" => raw_data_Sch9P_glutamine_L,
"Sch9_gtr1Delta" => raw_data_Sch9_gtr1Delta, "Sch9_glucose_starve" => raw_data_Sch9_glucose_starve,
"Sch9_glucose_relief" => raw_data_Sch9_glucose_relief, "Mig1_glucose_relief" => raw_data_Mig1_glucose_relief,
"Rib_rap" => raw_data_Rib_rap)

u_dict = Dict("sch9Delta_cAMP" => "cAMP", "Sch9P_glutamine_H" => "Sch9", "Sch9P_glutamine_L" => "Sch9",
"Sch9_gtr1Delta" => "Sch9","Sch9_glucose_starve" => "Sch9", "Sch9_glucose_relief" => "Sch9",
"Mig1_glucose_relief" => "Mig1", "Rib_rap" => "Rib")

pre_shift_dict = Dict("Snf1" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"cAMP" => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
"sch9Delta_cAMP" => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 0.0),
"Sch9P_glutamine_H" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0),
"Sch9P_glutamine_L" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0),
"Sch9_gtr1Delta" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0),
"Sch9_glucose_starve" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"Sch9_glucose_relief" => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
"Mig1_glucose_relief" => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
"Rib_rap" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))

post_shift_dict = Dict("Snf1" => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
"cAMP" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"sch9Delta_cAMP" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"Sch9P_glutamine_H" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"Sch9P_glutamine_L" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3),
"Sch9_gtr1Delta" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"Sch9_glucose_starve" => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
"Sch9_glucose_relief" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"Mig1_glucose_relief" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
"Rib_rap" => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))

pre_p_const_changes = Dict("sch9Delta_cAMP" => Dict("Sch9_T" => (Sch9_T => 0.0)),
"Sch9_gtr1Delta" => Dict("EGO_T" => (EGO_T => 0.0)))

post_p_const_changes = Dict("sch9Delta_cAMP" => Dict("Sch9_T" => (Sch9_T => 0.0)),
"Rib_rap" => Dict("TORC1_T" => (TORC1_T => 0.0)))

pre_p_var_changes = Dict("Sch9_gtr1Delta" => Dict("w_torc_ego" => 0.0,
"w_torc_egoin" => 0.0))