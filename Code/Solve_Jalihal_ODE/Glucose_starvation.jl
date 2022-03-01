# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")
# Pre-shift => Carbon, ATP = 1
p[findfirst(x->x=="Carbon", p_lookup_table)] = Carbon => 1.0
p[findfirst(x->x=="ATP", p_lookup_table)] = ATP => 1.0

include("ODE_methods.jl")
u0_SS = Steady_state_solver(p) # Returnerar steady state för parametrarna p

# Post-shift, Glucose addition => Carbon, ATP = 0
p[findfirst(x->x=="Carbon", p_lookup_table)] = Carbon => 0.0
p[findfirst(x->x=="ATP", p_lookup_table)] = ATP => 0.0

tspan = (0.0, 60.0) # [min]
sol = ODE_solver(u0_SS, tspan, p)

plot(sol, vars=Snf1, legend=false)
xlabel!("t [min]")
ylabel!("Snf1")
ylims!((0.0, 1.0))
title!("Glucose starvation, Carbon = 0, ATP = 0")