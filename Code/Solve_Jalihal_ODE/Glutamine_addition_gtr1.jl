# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")

# Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
p[findfirst(x->x=="EGO_T", p_lookup_table)] = EGO_T => 0.0
p[findfirst(x->x=="w_torc_ego", p_lookup_table)] = w_torc_ego => 0.0
p[findfirst(x->x=="w_torc_egoin", p_lookup_table)] = w_torc_egoin => 0.0

# Pre-shift => Glutamine_ext = 0
p[findfirst(x->x=="Glutamine_ext", p_lookup_table)] = Glutamine_ext => 0.0

include("ODE_methods.jl")
u0_SS = Steady_state_solver(p) # Returnerar steady state för parametrarna p

# High glutamine => Glutamine_ext = 1
p[findfirst(x->x=="Glutamine_ext", p_lookup_table)] = Glutamine_ext => 1.0

tspan = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, tspan, p)

plot(sol, vars=Sch9, legend=false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.0))
title!("Glutamine addition gtr1\\Delta, Glutamine_{ext} = 1.0")