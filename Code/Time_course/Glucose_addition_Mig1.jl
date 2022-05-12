# Creates the plot Mig1 under glucose addition with pre- and postshift variables

using DifferentialEquations
using ModelingToolkit
using Plots

default(titlefontsize=13)
default(dpi=300)

include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0 which is glucose starvation
include("../Model/ODE_methods.jl")

# Return the steady state for the parameters with glucose satarvation
u0_SS = Steady_state_solver(p_const, p_var, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) 


# Post-shift => ATP, Carbon = 1 which is clucos addition after sarvation
tspan = (0.0, 20.0) # [min]

# Solve the ODEs for the u0 from glucose addition over tspan
sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var) 

# Includes the experimental data for Mig1 glucose relief and plots the data as a scatter-plot
include("../../Data/exp_data_norm.jl")
scatter(
    t_Mig1_glucose_relief, data_Mig1_glucose_relief, markersize = 4.5, 
    markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8
    )

# Plot the results for Mig1
plot1 = plot!(sol, vars=log(10, Mig1/(1e-5 + 1.0-Mig1)), color = color_jalihal, lw=2.0, legend=false) 
xlabel!("t [min]")
ylabel!("log(nMig1/cMig1)")
ylims!((1.0, 1.6))
xticks!([0.0, 10.0, 20.0])
yticks!([1.0, 1.2, 1.4, 1.6])
title!("Glukostillsättning")
display(plot1)

file_name = "Mig1_gluc_add_reconstr"
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 
