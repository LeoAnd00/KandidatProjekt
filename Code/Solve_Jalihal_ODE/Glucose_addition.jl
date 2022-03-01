# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 

include("parameter_values.jl")
# Pre-shift => ATP, Carbon = 0
p[findfirst(x->x=="ATP", p_lookup_table)] = ATP => 0.0 
p[findfirst(x->x=="Carbon", p_lookup_table)] = Carbon => 0.0

include("ODE_methods.jl")
u0_SS = Steady_state_solver(p) # Returnerar steady state för parametrarna p

# Post-shift, Glucose addition => ATP, Carbon = 1
p[findfirst(x->x=="ATP", p_lookup_table)] = ATP => 1.0 
p[findfirst(x->x=="Carbon", p_lookup_table)] = Carbon => 1.0

tspan_cAMP = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, tspan_cAMP, p)

include("exp_data.jl")
scatter(t_wt_cAMP,data_wt_cAMP)

plot1 = plot!(sol, vars = cAMP, legend = false)
xlabel!("t [min]")
ylabel!("cAMP")
title!("Glucose addition, Carbon = 1.0")
display(plot1)

tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, tspan_Sch9, p)

include("exp_data.jl")
scatter(t_Sch9_glucose_relief,data_Sch9_glucose_relief)

plot2 = plot!(sol, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
title!("Glucose addition, Carbon = 1.0")
display(plot2)