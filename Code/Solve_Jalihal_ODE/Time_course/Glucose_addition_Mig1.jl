# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")

# Pre-shift => ATP, Carbon = 0
include("ODE_methods.jl")
u0_SS = Steady_state_solver(p, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift => ATP, Carbon = 1
tspan = (0.0, 20.0) # [min]
sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p)

include("exp_data.jl")
scatter(t_Mig1_glucose_relief, data__Mig1_glucose_relief, markersize = 5.5, markercolor="blue", markerstrokewidth=1.5)

p1 = plot!(sol, vars=log(10, Mig1/(1.0-Mig1)), color = "red", lw=2.5, legend=false)
xlabel!("t [min]")
ylabel!("log(Mig1/(Mig1_T-Mig1)")
ylims!((1.0, 1.6))
# yticks!([1.2, 1.4])
title!("Glucose addition")
display(p1)
