# 

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

include("../Model/ODE_methods.jl")
include("../../Data/exp_data_norm.jl")
include("../../Results/Parameter_values/Results_from_optimization.jl")
u0_SS = Steady_state_solver(p_const, p_var_alt_1, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift: Rapamycin treatment => TORC1_T = 0.0 
p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0

# Steady-state vid svältning utgör initialvärden för Glutamine addition
tspan = (0.0, 90.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var_alt_1)

plot1 = scatter(t_Rib_rap, data_Rib_rap, labels = "Data")

plot!(sol, vars=Rib/(1e-3+u0_SS[25]), color = Color_alt_1, lw=2.5, labels= label_alt_1) # 1e-2 ensures denom != 0

p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 1.0
u0_SS = Steady_state_solver(p_const, p_var_alt_2, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0
# Steady-state vid svältning utgör initialvärden för Glutamine addition
tspan = (0.0, 90.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var_alt_2)

plot!(plot1, sol, vars=Rib/(1e-3+u0_SS[25]), color = Color_alt_2, lw=2.5, labels= label_alt_2)

p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 1.0
u0_SS = Steady_state_solver(p_const, p_var_jalihal, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0
# Steady-state vid svältning utgör initialvärden för Glutamine addition
tspan = (0.0, 90.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var_jalihal)

plot!(plot1, sol, vars=Rib/(1e-3+u0_SS[25]), color = Color_jalihal, lw=2.5, labels= label_jalihal)

xlabel!("t [min]")
ylabel!("Rel. RPL32 mRNA") # Motsvarar konc. Rib relativt steady state
title!("Rapamycin treatment, TORC1_T = 0")
display(plot1)
savefig(plot1, pwd()*"/Results/Plots_for_new_and_old_paramter_values/Rapamycin.png")