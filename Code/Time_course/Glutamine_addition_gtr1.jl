# Creates plot for  Sch9 with a mutation gtr1_Delta under Glutamine addition

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")

# Add the mutant gtr1_Delta by setting EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 
p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0

# Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1, glutamine starvation
include("../Model/ODE_methods.jl")
# Returns the u0 (startvalues) from the pre-shift 
u0_SS = Steady_state_solver(p_const, p_var, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

# High glutamine => Glutamine_ext = 1, glutamine addition
tspan = (0.0, 30.0) # [min]
# Solve the ODEs for the post-shift changes
sol = ODE_solver(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Plot the solution with experimental data from a seperate file
plot1 = scatter(t_Sch9_gtr1Delta, data_Sch9_gtr1Delta)

plot!(sol, vars=Sch9, legend=false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
title!("Glutamine addition gtr1\\Delta, Glutamine_{ext} = 1.0")
display(plot1)