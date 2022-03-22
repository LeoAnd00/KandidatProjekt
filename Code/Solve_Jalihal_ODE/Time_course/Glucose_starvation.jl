# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")

# Pre-shift => Carbon, ATP = 1
include("ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift, Glucose starvation => Carbon, ATP = 0
tspan_Snf1 = (0.0, 61.0) # [min]
sol_Snf1 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Snf1, p_const, p_var)

include("exp_data.jl")

data = Minmaxnorm(data_Snf1)
plot1 = scatter(t_Snf1, data)

plot!(sol_Snf1, vars=Snf1, legend=false)
xlabel!("t [min]")
ylabel!("Snf1")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Snf1)*1.02))
title!("Glucose starvation")
display(plot1)

tspan_Sch9 = (0.0, 30.0) # [min]
sol_sch9 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

include("exp_data.jl")

data = Minmaxnorm(data_Sch9_glucose_starve)
plot2 = scatter(t_Sch9_glucose_starve, data)

plot!(sol_sch9, vars=Sch9, legend=false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Sch9)*1.02))
title!("Glucose starvation")
display(plot2)