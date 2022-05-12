# Creates the plot for Rapamycin treatment where Ribi is used to make the plot

using DifferentialEquations
using ModelingToolkit
using Plots

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")
# Pre-shift => No change in parameter values

include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) 

# Post-shift: Rapamycin treatment, which means TORC1_T = 0.0 
p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0

# Steady-state at starvation which is the initial values for Glutamine addition
tspan = (0.0, 90.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var) 

include("../../Data/exp_data_norm.jl")

plot10 = scatter(t_Rib_rap, data_Rib_rap, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

# Rel. RPL32 mRNA corresponds to (Rib divided by the pre-sfhit steady state value for Rib)
plot!(sol, vars=Rib/(1e-3+u0_SS[Get_index(u_lookup_table, "Rib(t)")]), legend=false, color=color_jalihal, lw=2.0) 
xlabel!("t [min]")
ylabel!("Rel. RPL32 mRNA") # Corresponds to Rib relativt steady state
title!("Rapamycinbehandling (TORC1 = 0)")
display(plot10)

file_name = "Rib_rap_treat_reconstr"
savefig(plot10, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot10, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 