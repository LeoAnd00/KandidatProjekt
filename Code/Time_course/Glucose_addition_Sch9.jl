# Creates the plot cAMP under glucose addition with pre- and postshift variables

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")
# Mutant Sch9_Delta => Sch9_T = 0, Creates a mutant by deleting Sch9_T
p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0

# Pre-shift => ATP, Carbon = 0, which is glucos starvation
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Return steady state parameters for glucos starvation

# Post-shift => ATP, Carbon = 1, which is glucos addition 
tspan = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)  # Solves the ODEs for the u0 from glucos addition over tspan

# Makes a scatter-plot from the experimental data from cAMP 
include("../../Data/exp_data_norm.jl")
scatter(t_sch9Delta_cAMP, data_sch9Delta_cAMP)

# Plot results from simulation and experimental data
plot1 = plot!(sol, vars=cAMP, legend=false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0, 1.1))
xlims!((-0.03, last(tspan)*1.01))
title!("Glucose addition Sch9\\Delta")
display(plot1)