# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")

# Pre-shift => Glutamine_ext = 0
include("ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady state för parametrarna p

# Low glutamine => Glutamine_ext = 0.3
tspan = (0.0, 30.0) # [min]
sol_low = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, p_var)

include("exp_data.jl")

data = Minmaxnorm(data_Sch9P_glutamine_L)
plot1 = scatter(t_Sch9P_glutamine_L, data)

plot!(sol_low, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
title!("Glutamine addition, Low Glutamine (0.3)")
display(plot1)

# High Glutamine => Glutamine_ext = 1.0
sol_high = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

include("exp_data.jl")

data = Minmaxnorm(data_Sch9P_glutamine_H)
plot2 = scatter(t_Sch9P_glutamine_H, data)

plot!(sol_high, vars = Sch9, label = "High Glutamine (1.0)", legend = false)
xlabel!("t [min]")
xlims!((-0.5, last(tspan)*1.02))
ylabel!("Sch9")
title!("Glutamine addition, High glutamine")
display(plot2)