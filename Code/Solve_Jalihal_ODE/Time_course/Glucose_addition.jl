# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")

# Pre-shift => ATP, Carbon = 0
include("ODE_methods.jl")
u0_SS = Steady_state_solver(p, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift, Glucose addition => ATP, Carbon = 1
tspan_cAMP = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_cAMP, p)

include("exp_data.jl")
scatter(t_wt_cAMP, data_wt_cAMP)

plot1 = plot!(sol, vars = cAMP, legend = false)
xlabel!("t [min]")
ylabel!("cAMP")
title!("Glucose addition")
display(plot1)

tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_Sch9, p)

scatter(t_Sch9_glucose_relief,data_Sch9_glucose_relief)

plot2 = plot!(sol, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
title!("Glucose addition")
display(plot2)