# Creates plot for high and low glutamine addition after starvation with Sch9 as a readout

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")

# Pre-shift => Glutamine_ext = 0, glutamine starvation
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returns steady state values for p from pre-shift

# Low glutamine => Glutamine_ext = 0.3
tspan = (0.0, 30.0) # [min]
# Solve the ODEs with low additon of glutamin as post shift
sol_low = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Creates the plot for low glutamine addition after starvation for Sch9
plot1 = scatter(t_Sch9P_glutamine_L, data_Sch9P_glutamine_L, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol_low, vars = Sch9, legend = false, color=color_jalihal)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
title!("Glutamine addition, Low Glutamine (0.3)")
display(plot1)

file_name = "Sch9_nit_add_low_reconstr"
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 

# High Glutamine => Glutamine_ext = 1.0
# Solves the same ODEs with high glutamin additon 
sol_high = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

# Plot the solution for high glutamine addition, with experimental data form seperate file
plot2 = scatter(t_Sch9P_glutamine_H, data_Sch9P_glutamine_H, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol_high, vars = Sch9, label = "High Glutamine (1.0)", legend = false, color=color_jalihal)
xlabel!("t [min]")
xlims!((-0.5, last(tspan)*1.02))
ylabel!("Sch9")
title!("Glutamine addition, High glutamine")
display(plot2)

file_name = "Sch9_nit_add_high_reconstr"
savefig(plot2, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot2, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 