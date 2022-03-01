# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")
# Pre-shift => Glutamine_ext = 0
p[findfirst(x->x=="Glutamine_ext", p_lookup_table)] = Glutamine_ext => 0.0

include("ODE_methods.jl")
u0_SS = Steady_state_solver(p) # Returnerar steady state för parametrarna p

# Low glutamine => Glutamine_ext = 0.3
p[findfirst(x->x=="Glutamine_ext", p_lookup_table)] = Glutamine_ext => 0.3 

tspan = (0.0, 30.0) # [min]

sol_low = ODE_solver(u0_SS, tspan, p)

# High Glutamine => Glutamine_ext = 1.0
p[findfirst(x->x=="Glutamine_ext", p_lookup_table)] = Glutamine_ext => 1.0

sol_high = ODE_solver(u0_SS, tspan, p)

plot(sol_low, vars = Sch9, label = "Low Glutamine (0.3)")
plot!(sol_high, vars = Sch9, label = "High Glutamine (1.0)", legend=:top)
xlabel!("t [min]")
ylabel!("Sch9")
title!("Glutamine addition, Glutamine_ext = 0.3")