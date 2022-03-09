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
u0_SS_sit4 = Steady_state_solver(p, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) 

# Post-shift, Rapamycin treatment Sit4 => TORC1_T = 0 => no change in parameter values
tspan_sit4 = (0.0, 150.0) # [min]
sol_sit4 = ODE_solver(u0_SS_sit4, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_sit4, p_sit4)

########## wt #########
include("parameter_values.jl")
p_wt = p

# Pre-shift wt => no change in parameter values

u0_SS_wt = Steady_state_solver(p_wt, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift, Rapamycin treatment wt => TORC1_T = 0
p_wt[findfirst(x->x=="TORC1_T", p_lookup_table)] = TORC1_T => 0.0

tspan_wt = (90.0, 150.0)
sol_wt = ODE_solver(u0_SS_wt, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_wt, p_wt)

plot(sol_sit4, vars = 1-Gln3, label = "Sit4 (1-Gln3)", lw = 4, ls=:dash, legend=:right)
plot!(sol_wt, vars = Gln3, label = "wt (Gln3)", lw = 4, ls=:dash, color = "red")
p1 = plot!([0.0, 90.0], [u0_SS_wt[16], u0_SS_wt[16]], ls=:dash, lw=4, color = "red", label="")

xlabel!("t [min]")
ylabel!("Gln3") 
xlims!((0, 150))
ylims!((0, 1.0))
title!("Rapamycin treatment, TORC1_T = 0")
display(p1)