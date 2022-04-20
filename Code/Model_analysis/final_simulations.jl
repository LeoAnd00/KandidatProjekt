include("../Model/ODE_methods.jl")
include("../Model/parameter_values.jl")

using Plots

blue_Jali = RGB(50/255,150/255,245/255)
red_alt1 = RGB(250/250,140/255,0/255)
orange_alt2 = RGB(150/250,200/250,48/250)

default(dpi = 300)

###### Glukossvältning ######
u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))
sol_gluc_starve =  ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_starvation]), (0.0, 30.0), p_conc)

plot1 = plot(sol_gluc_starve, vars = [Gln3, Mig1], title = "Glukossvältning", legend=:right) # Gln3 & Cyr1 beror på Snf1
ylims!((0.0, 1.0))
display(plot1)

###### Kvävesvältning ######
u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_nitrogen_starvation]))
sol_nit_starve =  ODE_solver(u0_SS, last.(nutrient_shifts[index_nitrogen_starvation]), (0.0, 10.0), p_conc)

plot_nit_starve = plot(sol_nit_starve, vars = [Snf1], title = "Kvävesvältning")
ylims!((0.0, 1.0))
display(plot_nit_starve)

###### Glukostillsättning ######
u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_addition]))
sol_gluc_add =  ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_addition]), (0.0, 20.0), p_conc)

plot_gluc_add = plot(sol_gluc_add, vars = [cAMP, PKA, Sak], title = "Glukostillsättning", legend=:right, color=[blue_Jali red_alt1 orange_alt2]) # Oscillerande ämnen
display(plot_gluc_add)

###### början av glukossvältning (låg) ######
u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))
sol =  ODE_solver(u0_SS, (Carbon => 0.3, ATP => 0.3, Glutamine_ext => 1.0), (0.0, 10.0), p_conc)

###### Kvävetillsättning ######
u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_high_glutamine]))
sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_high_glutamine]), (0.0, 34.0), p_conc)

###### Snf1 deletion ######

include("../Model/parameter_values.jl")

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))
sol_gluc_starve =  ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_starvation]), (0.0, 10.0), p_conc)

p_conc[Get_index(p_conc_lookup_table, "Snf1_T")] = Snf1_T => 0.0

u0_SS = Steady_state_solver(p_conc, first.(nutrient_shifts[index_glucose_starvation]))
sol = ODE_solver(u0_SS, last.(nutrient_shifts[index_glucose_starvation]), (0.0, 10.0), p_conc)

plot3 = plot(sol, vars = [Cyr1], title = "Snf1\$\\Delta\$")
plot!(sol_gluc_starve, vars = [Cyr1], label = "Cyr1_wt")
display(plot3)