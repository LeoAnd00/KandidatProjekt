# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")
# Mutant Sch9_Delta => Sch9_T = 0
p[findfirst(x->x=="Sch9_T", p_lookup_table)] = Sch9_T => 0.0 

# Pre-shift => ATP, Carbon = 0
include("ODE_methods.jl")
u0_SS = Steady_state_solver(p, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift => ATP, Carbon = 1
tspan = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p)

include("exp_data.jl")
scatter(t_sch9Delta_cAMP, data_sch9Delta_cAMP)

p1 = plot!(sol, vars=cAMP, legend=false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0, 1.0))
xlims!((0, 3.05))
title!("Glucose addition Sch9\\Delta")
display(p1)