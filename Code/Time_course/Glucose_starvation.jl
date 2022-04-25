# Creates the plots for Glucose starvation for Snf1 and Sch9 separtly sincethey have diffrent time stamps 

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Pre-shift => Carbon, ATP = 1, glucose abundence 
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returns the steady state value for p

# Post-shift, Glucose starvation => Carbon, ATP = 0
tspan_Snf1 = (0.0, 61.0) # [min]
# Solves the ODEs with parametrs for Glucose starvation for Snf1
sol_Snf1 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Snf1, p_const, p_var) 

include("../../Data/exp_data_norm.jl")

# Creates the plot for Snf1 with experimental data from another file
plot1 = scatter(t_Snf1, data_Snf1)
plot!(sol_Snf1, vars=Snf1, legend=false)
xlabel!("t [min]")
ylabel!("Snf1")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Snf1)*1.02))
title!("Glucose starvation")
display(plot1)

tspan_Sch9 = (0.0, 30.0) # [min]
# Solves the ODEs for Sch9 in the same way as Snf1 but for a diffrent tspan
sol_sch9 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

# Plot solution with ex data from a diffrent file
plot2 = scatter(t_Sch9_glucose_starve, data_Sch9_glucose_starve)
plot!(sol_sch9, vars=Sch9, legend=false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Sch9)*1.02))
title!("Glucose starvation")
display(plot2)