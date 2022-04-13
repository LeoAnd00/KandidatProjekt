# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")
# Mutant Sch9_Delta => Sch9_T = 0
p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0

# Pre-shift => ATP, Carbon = 0
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift => ATP, Carbon = 1
tspan = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")

scatter(t_sch9Delta_cAMP, data_sch9Delta_cAMP)

plot1 = plot!(sol, vars=cAMP, legend=false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0, 1.1))
xlims!((-0.03, last(tspan)*1.01))
title!("Glucose addition Sch9\\Delta")
display(plot1)