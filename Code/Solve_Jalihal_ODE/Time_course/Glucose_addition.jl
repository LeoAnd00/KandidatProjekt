# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")

# Pre-shift => ATP, Carbon = 0
include("ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift, Glucose addition => ATP, Carbon = 1
tspan_cAMP = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_cAMP, p_const, p_var)

include("exp_data.jl")

data = Minmaxnorm(data_cAMP)
plot1 = scatter(t_cAMP, data)

plot!(sol, vars = cAMP, legend = false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0.0, 1.02))
xlims!((-0.02, last(tspan_cAMP)*1.02))
title!("Glucose addition")
display(plot1)

tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

data = Minmaxnorm(data_Sch9_glucose_relief)
plot2 = scatter(t_Sch9_glucose_relief, data)

plot!(sol, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.02))
xlims!((-0.3, last(tspan_Sch9)*1.02))
title!("Glucose addition")
display(plot2)
