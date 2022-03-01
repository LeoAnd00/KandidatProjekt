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

include("exp_data.jl")
scatter(t_Sch9P_glutamine_L,data_Sch9P_glutamine_L)

plot1 = plot!(sol_low, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
title!("Glutamine addition, Low Glutamine (0.3)")
display(plot1)

# High Glutamine => Glutamine_ext = 1.0
p[findfirst(x->x=="Glutamine_ext", p_lookup_table)] = Glutamine_ext => 1.0

sol_high = ODE_solver(u0_SS, tspan, p)

include("exp_data.jl")
scatter(t_Sch9P_glutamine_H,data_Sch9P_glutamine_H)

plot2 = plot!(sol_high, vars = Sch9, label = "High Glutamine (1.0)",legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
title!("Glutamine addition, High glutamine")
display(plot2)