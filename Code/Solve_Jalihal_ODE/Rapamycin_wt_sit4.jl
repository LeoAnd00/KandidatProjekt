# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

########## Sit4 #########
include("parameter_values.jl")
# Deklarerar separat parametervektor för sit4
p_sit4 = p 

# Pre-shift Sit4 => TORC1_T = 0
p_sit4[findfirst(x->x=="TORC1_T", p_lookup_table)] = TORC1_T => 0.0

include("ODE_methods.jl")
u0_SS_sit4 = Steady_state_solver(p) # Returnerar steady state för parametrarna p

# Post-shift, Rapamycin treatment Sit4 => TORC1_T = 0 => no change in parameter values
tspan = (0.0, 90.0) # [min]

sol_sit4 = ODE_solver(u0_SS_sit4, tspan, p_sit4)

########## wt #########
include("parameter_values.jl")
p_wt = p

# Pre-shift wt => no change in parameter values

u0_SS_wt = Steady_state_solver(p_wt) # Returnerar steady state för parametrarna p

# Post-shift, Rapamycin treatment wt => TORC1_T = 0
p_wt[findfirst(x->x=="TORC1_T", p_lookup_table)] = TORC1_T => 0.0

sol_wt = ODE_solver(u0_SS_wt, tspan, p_wt)

plot(sol_sit4, vars = 1-Gln3, label = "Sit4 (1-Gln3)", linewidth = 4, linestyle=:dash, legend=:right)
plot!(sol_wt, vars = Gln3, label = "wt (Gln3)", linewidth = 4, linestyle=:dash)
xlabel!("t [min]")
ylabel!("Gln3") 
ylims!((0, 1.0))
title!("Rapamycin treatment, TORC1_T = 0")