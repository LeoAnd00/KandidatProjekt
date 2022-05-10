include("simulation_methods.jl")

###### Glucose starvation ######
sol_array_gluc_starve = Generate_sol_array(index_glucose_starvation, (0.0, 40.0)) # Generates array of simulation

plot_simulation_result(true, sol_array_gluc_starve, Mig1, (0.0, 20.0), (0.8, 1.0), "Glukossvältning", "Mig1_gluc_starve", :topright) 
plot_simulation_result(true, sol_array_gluc_starve, Snf1, (0.0, 1.0), "Glukossvältning", "Snf1_gluc_starve", :right) 
plot_simulation_result(true, sol_array_gluc_starve, Gis1, (0.0, 10.0), (0.0, 1.0), "Glukossvältning", "Gis1_gluc_starve", :topright) 
plot_simulation_result(true, sol_array_gluc_starve, Gln3, (0.0, 1.0), "Glukossvältning", "Gln3_gluc_starve", :right) 
plot_simulation_result(true, sol_array_gluc_starve, Dot6, (0.0, 4.0), (0.8, 1.0), "Glukossvältning", "Dot6_gluc_starve", :topright) 


###### Nitrogen starvation ######
sol_array_nit_starve = Generate_sol_array(index_nitrogen_starvation, (0.0, 40.0))

plot_simulation_result(true, sol_array_nit_starve, Snf1, (0.0, 15.0), (0.0, 0.2), "Kvävesvältning", "Snf1_nit_starve", :right) 
plot_simulation_result(true, sol_array_nit_starve, TORC1, (0.0, 1.0), "Kvävesvältning", "TORC1_nit_starve", :right) 
plot_simulation_result(true, sol_array_nit_starve, Sch9, (0.0, 1.0), "Kvävesvältning", "Sch9_nit_starve", :right) 
plot_simulation_result(true, sol_array_nit_starve, Gln3, (0.0, 1.0), "Kvävesvältning", "Gln3_nit_starve" ,:right) 
plot_simulation_result(true, sol_array_nit_starve, Rtg13, (0.0, 1.0), "Kvävesvältning", "Rtg13_nit_starve" ,:bottomright) 
plot_simulation_result(true, sol_array_nit_starve, Dot6, (0.0, 1.0), "Kvävesvältning", "Dot6_nit_starve" ,:right) 


###### Glucose addition ######
sol_array_gluc_add = Generate_sol_array(index_glucose_addition, (0.0, 15.0))

plot_simulation_result(true, sol_array_gluc_add, cAMP, (0.0, 15.0), (0.0, 1.05), "Glukostillsättning", "cAMP_gluc_add", :right) 
plot_simulation_result(true, sol_array_gluc_add, PKA, (0.0, 15.0), (0.0, 1.05), "Glukostillsättning", "PKA_gluc_add", :bottomright) 
plot_simulation_result(true, sol_array_gluc_add, Snf1, (0.0, 10.0), (0.0, 1.0), "Glukostillsättning", "Snf1_gluc_add", :right) 
plot_simulation_result(true, sol_array_gluc_add, Gis1, (0.0, 7.0), (0.0, 1.0), "Glukostillsättning", "Gis1_gluc_add", :right) 

###### Gtr1Delta | EGO deletion, Nitrogen addition ######
sol_array_nit_add = Generate_sol_array(index_high_glutamine, (0.0, 40.0))
sol_array_gtr1Delta = Generate_sol_array(index_high_glutamine, (0.0, 40.0), [EGO_T => 0.0])
plot_wt_mut_results(true, sol_array_nit_add, sol_array_gtr1Delta, Gln3, (0.0, 25.0), (0.0, 1.0), "Kvävetillsättning, gtr1\$\\Delta\$", "Gln3_gtr1Delta_nit_add", :right)

###### Snf1 deletion, Glucose starvation ######
sol_array_Snf1Delta = Generate_sol_array(index_glucose_starvation, (0.0, 1.0), [Snf1_T => 0.0])

plot_wt_mut_results(true, sol_array_gluc_starve, sol_array_Snf1Delta, Gln3, (0.0, 0.5), (0.0, 0.8), "Glukossvältning, Snf1\$\\Delta\$", "Cyr1_gluc_add", :right)

plot_wt_mut_results(true, sol_array_gluc_starve, sol_array_Snf1Delta, Cyr1, (0.0, 0.5), (0.0, 0.8), "Glukossvältning, Snf1\$\\Delta\$", "Cyr1_gluc_add", :right)