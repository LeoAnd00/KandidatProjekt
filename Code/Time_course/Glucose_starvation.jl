# Creates the plots for Glucose starvation for Snf1 and Sch9 separtaley since they have different time spans

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Pre-shift => Carbon, ATP = 1, glucose abundence 
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returns the steady state value for p

# Post-shift, Glucose starvation => Carbon, ATP = 0
tspan_Snf1 = (0.0, 61.0) # [min]

# Solves the ODEs under glucose starvation
sol_Snf1 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Snf1, p_const, p_var) 

include("../../Data/exp_data_norm.jl")

# Creates the plot for Snf1 with experimental data from another file
plot5 = scatter(t_Snf1, data_Snf1, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol_Snf1, vars=Snf1, legend=false, color=color_jalihal, lw=2.0)
xlabel!("t [min]")
ylabel!("Snf1")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Snf1)*1.02))
title!("Glukossvältning")
display(plot5)

file_name = "Snf1_gluc_starve_reconstr"
savefig(plot5, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot5, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 

tspan_Sch9 = (0.0, 30.0) # [min]

# Solves the ODEs for Sch9 in the same way as for Snf1 but with a different tspan
sol_sch9 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

# Plot solution with exp data from a different file
plot6 = scatter(t_Sch9_glucose_starve, data_Sch9_glucose_starve, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)
plot!(sol_sch9, vars=Sch9, legend=false, color=color_jalihal)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan_Sch9)*1.02))
title!("Glukossvältning")
display(plot6)

file_name = "Sch9_gluc_starve_reconstr"
savefig(plot6, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot6, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf")