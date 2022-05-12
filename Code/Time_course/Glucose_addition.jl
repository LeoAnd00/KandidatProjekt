# Creates the plot for cAMP and Sch9 for the wildtype

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0, glucose starvation
include("../Model/ODE_methods.jl")

# Returns steady state parameters that are used as u0 (start values for solving the ODE)
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 

# Post-shift, Glucose addition => ATP, Carbon = 1, glucose additon after starvation
tspan_cAMP = (0.0, 3.0) # [min]

# Solve the ODEs under the give time span and values on p
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_cAMP, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Create plot for cAMP and Sch9 separately with experimental data for each
plot3 = scatter(t_cAMP, data_cAMP, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol, vars = cAMP, legend = false, color=color_jalihal, lw=2.0, show=false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0.0, 1.02))
xlims!((-0.02, last(tspan_cAMP)*1.02))
title!("Glukostillsättning")
display(plot3)

file_name = "cAMP_gluc_add_reconstr"
savefig(plot3, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot3, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 

# Separate solution for Sch9 since it needs a bigger time span
tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

# Plot the solutions for Sch9
plot4 = scatter(
    t_Sch9_glucose_relief, data_Sch9_glucose_relief, markersize = 4.5, 
    markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8
    )
plot!(sol, vars = Sch9, legend = false, color=color_jalihal, lw=2.0)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.02))
xlims!((-0.3, last(tspan_Sch9)*1.02))
title!("Glukostillsättning")
display(plot4)

file_name = "Snf1_gluc_starve_reconstr"
savefig(plot4, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot4, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 
