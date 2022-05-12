include("Model_analysis/simulation_methods.jl")
sol_array_gluc_add = Generate_sol_array(index_glucose_addition, (0.0, 15.0))
plot_simulation_result(
    true, sol_array_gluc_add, cAMP, (0.0, 15.0), (0.0, 1.05), "Glukostillsättning", "cAMP_gluc_add", :right
    ) 
plot_simulation_result(
    true, sol_array_gluc_add, PKA, (0.0, 15.0), (0.0, 1.05), "Glukostillsättning", "PKA_gluc_add", :bottomright
    ) 
include("Time_course/plot_all_reconstr_in_subplot.jl")
include("Model_analysis/plot_all_sims_in_subplot.jl")