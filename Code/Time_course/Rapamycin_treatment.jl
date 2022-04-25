# Creates the plot for Rapamycin treatment where Ribi is used to make the plot

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")
# Pre-shift => No change in parameter values

include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returns steady state for the parameters p

# Post-shift: Rapamycin treatment which means TORC1_T = 0.0 
p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0

# Steady-state at starvation which is the initial values for Glutamine addition
tspan = (0.0, 90.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var) # Solve all ODEs

# Rel. RPL32 mRNA corresponds to (Rib/steady-state value for Rib)
include("../../Data/exp_data_norm.jl")

# Plot the results from the simulations and data from a seperate file
plot1 = scatter(t_Rib_rap, data_Rib_rap)

plot!(sol, vars=Rib/(1e-3+u0_SS[25]), legend=false) # 1e-2 ensures denom != 0
xlabel!("t [min]")
ylabel!("Rel. RPL32 mRNA") # Corresponds to Rib relativt steady state
title!("Rapamycin treatment, TORC1_T = 0")
display(plot1)