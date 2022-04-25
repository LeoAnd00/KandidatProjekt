# Creates the plot for cAMP and Sch9 withour mutants

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0, glucose starvation
include("../Model/ODE_methods.jl")

# Returns steady stat parameters that are used as u0 (start values for solving the ODE)
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 

# Post-shift, Glucose addition => ATP, Carbon = 1, glucose additon after starvation
tspan_cAMP = (0.0, 3.0) # [min]
# Solve the ODEs under the give time span and values on p
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_cAMP, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Creat plot for cAMP and Sch9 separtly with experimental data in each
plot1 = scatter(t_cAMP, data_cAMP)

plot!(sol, vars = cAMP, legend = false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0.0, 1.02))
xlims!((-0.02, last(tspan_cAMP)*1.02))
title!("Glucose addition")
display(plot1)

# Separat solever for Sch9 since it needs a bigger time stamp
tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

# Plot the solutions for Sch9
plot2 = scatter(t_Sch9_glucose_relief, data_Sch9_glucose_relief)
plot!(sol, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.02))
xlims!((-0.3, last(tspan_Sch9)*1.02))
title!("Glucose addition")
display(plot2)
