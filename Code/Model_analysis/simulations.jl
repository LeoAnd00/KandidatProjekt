include("../Model/ODE_methods.jl")
include("../Model/parameter_values.jl")

using Plots

###### Glukossvältning ######

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))

sol_gluc_starve =  ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_starvation]), (0.0, 10.0), p_conc)
plot1 = plot(sol_gluc_starve, vars = [Snf1, Gln3, Cyr1, Mig1], title = "Glukossvältning") # Gln3 & Cyr1 beror på Snf1
display(plot1)
plot1_2 = plot(sol_gluc_starve, vars = [Dot6, Gis1, Rib, Sch9], title = "Glukossvältning")
display(plot1_2)

###### Kvävesvältning ######

# Gln3, Dot6, Rtg13 bör defosforyleras, men gör ej det

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_nitrogen_starvation]))

sol_nit_starve =  ODE_solver(u0_SS, last.(nutrient_shifts[index_nitrogen_starvation]), (0.0, 50.0), p_conc)
plot2 = plot(sol_nit_starve, vars = [TORC1, Sch9, Gln3], title = "Kvävesvältning", legend=:right)
display(plot2)
plot2_2 = plot(sol_nit_starve, vars = [Rtg13, Dot6, Snf1], title = "Kvävesvältning")
display(plot2_2)

####### Glukostillsättning #######

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_addition]))

sol_gluc_add =  ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_addition]), (0.0, 10.0), p_conc)
# plot3 = plot(sol_gluc_add, vars = [Snf1, Gis1, Mig1, PKA, Sch9], title = "Glukostillsättning", legend=:right)
# plot3 = plot(sol_gluc_add, vars = [PKA, PDE, cAMP, Cyr1, Ras, Sak], title = "Glukostillsättning", legend=:right) # Oscillerande ämnen
plot3 = plot(sol_gluc_add, vars = [cAMP, PKA, Sak], title = "Glukostillsättning", legend=:right) # Oscillerande ämnen

display(plot3)

####### början av glukossvältning (låg) #######

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))

sol =  ODE_solver(u0_SS, (Carbon => 0.3, ATP => 0.3, Glutamine_ext => 1.0), (0.0, 10.0), p_conc)
plot_y = plot(sol, vars = [Gis1, PKA, Snf1], title = "början av glukossvältning (låg glukoshalt)", legend=:right) # Oscillerande ämnen
display(plot_y)

###### Kvävetillsättning ######

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_high_glutamine]))

sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_high_glutamine]), (0.0, 34.0), p_conc)
plot_x = plot(sol, vars = [TORC1, Gln3, Rtg13], legend=:right, title = "Kvävetillsättning")
display(plot_x)
plot_x2 = plot(sol, vars = [Dot6, Rib], legend=:right, title = "Kvävetillsättning")
display(plot_x2)

###### EGO deletion #######

p_conc[Get_index(p_conc_lookup_table, "EGO_T")] = EGO_T => 0.0

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_nitrogen_starvation]))

sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_nitrogen_starvation]), (0.0, 50.0), p_conc)

plot_T = plot(sol, vars = [Gln3], title = "EGO\$\\Delta\$", legend=:right)
plot!(sol_nit_starve, vars = [Gln3], label = "Gln3_wt")
display(plot_T)

###### Snf1 deletion ######

# Gln3 fosforyleras ej, vilket är förväntat.
# Cyr1 förväntas sluta defosforyleras, men den uppvisar samma beteende

include("../Model/parameter_values.jl")

p_conc[Get_index(p_conc_lookup_table, "Snf1_T")] = Snf1_T => 0.0

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))

sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_starvation]), (0.0, 10.0), p_conc)

plot5 = plot(sol, vars = [Gln3], title = "Snf1\$\\Delta\$")
plot!(sol_gluc_starve, vars = [Gln3], label = "Gln3_wt", legend=:right)
display(plot5)

plot5_2 = plot(sol, vars = [Cyr1], title = "Snf1\$\\Delta\$")
plot!(sol_gluc_starve, vars = [Cyr1], label = "Cyr1_wt")
display(plot5_2)

###### Sak1 deletion ######

# include("../Time_course/parameter_values.jl")

# p_conc[Get_index(p_conc_lookup_table, "Sak_T")] = Sak_T => 0.0

# u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))

# sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_starvation]), (0.0, 30.0), p_conc)

# plot4 = plot(sol, vars = [Snf1], title = "Sak\$\\Delta\$")
# plot!(sol_gluc_starve, vars = [Snf1], label = "wt")
# display(plot4)

###### Sit4 deletion ######

# include("../Time_course/parameter_values.jl")

# p_conc[Get_index(p_conc_lookup_table, "TORC1_T")] = TORC1_T => 0.0
# p_conc[Get_index(p_conc_lookup_table, "w_gln_sit")] = w_gln_sit => 0.0

# u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_nitrogen_starvation]))

# sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_nitrogen_starvation]), (0.0, 30.0), p_conc)

# plot6 = plot(sol, vars = [Gln3], title = "sit4\$\\Delta\$")
# plot!(sol_nit_starve, vars = [Gln3], label = "Gln3_wt")
# display(plot6)

# plot7 = plot(sol, vars = [Rtg13], title = "sit4\$\\Delta\$")
# plot!(sol_nit_starve, vars = [Rtg13], label = "Rtg13_wt")
# display(plot7)

# plot8 = plot(sol, vars = [Sch9], title = "sit4\$\\Delta\$")
# plot!(sol_nit_starve, vars = [Sch9], label = "Sch9_wt")
# display(plot8)