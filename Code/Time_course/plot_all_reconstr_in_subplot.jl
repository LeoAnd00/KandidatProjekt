# Plots all the reconstructed plots in a subplot 

using DifferentialEquations
using ModelingToolkit
using Plots

default(titlefontsize=5)
default(xtickfont=4)
default(ytickfont=4)
default(guidefont=6)
default(ylabelfontsize=4)
default(xlabelfontsize=5)
default(markersize=1.4)
default(lw=1.1)

##### Glucose addition Mig1 #####
include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0 which is glucose starvation
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0))


# Post-shift => ATP, Carbon = 1 which is clucos addition after sarvation
tspan = (0.0, 20.0) # [min]
sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

# Includes the experimental data for Mig1 glucos relief and presnet as a scatter-plot
include("../../Data/exp_data_norm.jl")
scatter(t_Mig1_glucose_relief, data_Mig1_glucose_relief, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

# Plot the results for Mig1
plot1 = plot!(sol, vars=log(10, Mig1/(1e-5 + 1.0-Mig1)), color = color_jalihal, legend=false) 
xlabel!("")
ylabel!("log(nMig1/cMig1)")
ylims!((1.0, 1.6))
xticks!([0.0, 10.0, 20.0])
yticks!([1.0, 1.2, 1.4, 1.6])
title1 = "Glukostillsättning"

##### Glucose addition Sch9_delta #####
include("../Model/parameter_values.jl")
# Mutant Sch9_Delta => Sch9_T = 0, Creates a mutant by deleting Sch9_T
p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0

# Pre-shift => ATP, Carbon = 0, which is glucos starvation
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 

# Post-shift => ATP, Carbon = 1, which is glucos addition 
tspan = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var) 

# Makes a scatter-plot from the experimental data from cAMP 
scatter(t_sch9Delta_cAMP, data_sch9Delta_cAMP, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

# Plot results from simulation and experimental data
plot2 = plot!(sol, vars=cAMP, legend=false, color=color_jalihal)
xlabel!("")
ylabel!("cAMP")
ylims!((0, 1.1))
xlims!((-0.03, last(tspan)*1.01))
yticks!([0.0, 0.50, 1.00])
title2="Glukostillsättning, Sch9\\Delta"

##### Glucose addition #####
include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0, glucose starvation
# Returns steady stat parameters that are used as u0 (start values for solving the ODE)
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 

# Post-shift, Glucose addition => ATP, Carbon = 1, glucose additon after starvation
tspan_cAMP = (0.0, 3.0) # [min]
# Solve the ODEs under the give time span and values on p
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_cAMP, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Creat plot for cAMP and Sch9 separately with experimental data in each
plot3 = scatter(t_cAMP, data_cAMP, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol, vars = cAMP, legend = false, color=color_jalihal, show=false)
xlabel!("")
ylabel!("cAMP")
ylims!((0.0, 1.1))
xlims!((-0.02, last(tspan_cAMP)*1.02))
yticks!([0.0, 0.50, 1.00])
title3="Glukostillsättning"

tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

# Plot the solutions for Sch9
plot4 = scatter(t_Sch9_glucose_relief, data_Sch9_glucose_relief,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol, vars = Sch9, legend = false, color=color_jalihal)
xlabel!("")
ylabel!("Sch9")
ylims!((0.0, 1.1))
yticks!([0.0, 0.50, 1.00])
xlims!((-0.3, last(tspan_Sch9)*1.02))
title4="Glukostillsättning"

##### Glucose starvation #####
# Pre-shift => Carbon, ATP = 1, glucose abundance 
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))

# Post-shift, Glucose starvation => Carbon, ATP = 0
tspan_Snf1 = (0.0, 61.0) # [min]
# Solves the ODEs with parametrs for Glucose starvation for Snf1
sol_Snf1 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Snf1, p_const, p_var) 

# Creates the plot for Snf1 with experimental data from another file
plot5 = scatter(t_Snf1, data_Snf1,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol_Snf1, vars=Snf1, legend=false, color=color_jalihal)
xlabel!("")
ylabel!("Snf1")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Snf1)*1.02))
yticks!([0.0, 0.50, 1.00])
xticks!([0.0, 20.0, 40.0, 60.0])
title5="Glukossvältning"

tspan_Sch9 = (0.0, 30.0) # [min]
# Solves the ODEs for Sch9 in the same way as Snf1 but for a diffrent tspan
sol_sch9 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

# Plot solution with ex data from a different file
plot6 = scatter(t_Sch9_glucose_starve, data_Sch9_glucose_starve,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol_sch9, vars=Sch9, legend=false, color=color_jalihal)
xlabel!("")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Sch9)*1.02))
yticks!([0.0, 0.50, 1.00])
title6="Glukossvältning"

##### Glutamine addition gtr1_Delta #####
p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0

# Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1, glutamine starvation
# Returns the u0 (startvalues) from the pre-shift 
u0_SS = Steady_state_solver(p_const, p_var, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

# High glutamine => Glutamine_ext = 1, glutamine addition
tspan = (0.0, 30.0) # [min]
# Solve the ODEs for the post-shift changes
sol = ODE_solver(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const, p_var)

# Plot the solution with experimental data from a seperate file
plot7 = scatter(t_Sch9_gtr1Delta, data_Sch9_gtr1Delta,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol, vars=Sch9, legend=false, color=color_jalihal)
xlabel!("")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
yticks!([0.0, 0.50, 1.00])
title7="Kvävetillsättning, gtr1\\Delta"

##### Glutamine addition #####
include("../Model/parameter_values.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) 

# Low glutamine => Glutamine_ext = 0.3
tspan = (0.0, 30.0) # [min]
# Solve the ODEs with low additon of glutamin as post shift
sol_low = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, p_var)

# Creates the plot for low glutamine addition after starvation for Sch9
plot8 = scatter(t_Sch9P_glutamine_L, data_Sch9P_glutamine_L,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol_low, vars = Sch9, legend = false, color=color_jalihal)
xlabel!("")
ylabel!("Sch9")
ylims!((0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
yticks!([0.0, 0.50, 1.00])
title8="Kvävetillsättning (låg, 0.3)"

sol_high = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

# Plot the solution for high glutamine addition, with experimental data form seperate file
plot9 = scatter(t_Sch9P_glutamine_H, data_Sch9P_glutamine_H,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol_high, vars = Sch9, legend = false, color=color_jalihal)
xlabel!("t [min]")
xlims!((-0.5, last(tspan)*1.02))
ylabel!("Sch9")
ylims!((0.0, 1.1))
yticks!([0.0, 0.50, 1.00])
title9="Kvävetillsättning (hög, 1.0)"

##### Rapamycin treatment #####
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) 

# Post-shift: Rapamycin treatment which means TORC1_T = 0.0 
p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0

# Steady-state at starvation which is the initial values for Glutamine addition
tspan = (0.0, 90.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var) 

# Rel. RPL32 mRNA corresponds to (Rib/steady-state value for Rib)

# Plot the results from the simulations and data from a seperate file
plot10 = scatter(t_Rib_rap, data_Rib_rap,  markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol, vars=Rib/(1e-3+u0_SS[Get_index(u_lookup_table, "Rib(t)")]), legend=false, color=color_jalihal, )
xlabel!("t [min]")
ylabel!("Rel. RPL32 mRNA") # Corresponds to Rib relative to it's steady-state value
yticks!([0.0, 0.50, 1.00])
ylims!((0.0, 1.1))
title10="Rapamycin (TORC1 = 0)"

title11="Vildtyp"
title12="Mutanter"

##### Perturbation #####
include("../Perturbation/perturbation_experiments.jl")

plot_list = [plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12]
l = @layout [
    grid(5,2){0.6w} grid(2,1)
]

titles=[title1, title2, title3, title4, title5, title6, title7, title8, title9, title10, title11, title12]

plot_all = plot(
    plot_list..., layout=l, legend = false, 
    title=["\\($i\\) "*title for j in 1:1, (i, title) in enumerate(titles)]
)
display(plot_all)

file_name = "reconstructed_subplot"
savefig(plot_all, pwd()*"/Results/"*file_name*".png")    
savefig(plot_all, pwd()*"/Results/"*file_name*".pdf") 