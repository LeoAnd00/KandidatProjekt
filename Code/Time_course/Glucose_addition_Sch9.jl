# Creates the plot cAMP under glucose addition with pre- and postshift variables

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Mutant Sch9_Delta => Sch9_T = 0, Creates a mutant by setting total amount of the protein to zero
p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0

# Pre-shift => ATP, Carbon = 0, which is glucose starvation
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 

# Post-shift => ATP, Carbon = 1, which is glucose addition 
tspan = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var) 

# Makes a scatter-plot from the experimental data from cAMP 
include("../../Data/exp_data_norm.jl")
scatter(t_sch9Delta_cAMP, data_sch9Delta_cAMP, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

# Plot results from simulation and experimental data
plot2 = plot!(sol, vars=cAMP, legend=false, color=color_jalihal, lw=2.0)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0, 1.1))
xlims!((-0.03, last(tspan)*1.01))
title!("Glukostillsättning, Sch9\\Delta")
display(plot2)

file_name = "cAMP_Sch9delta_gluc_starve_reconstr"
savefig(plot2, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot2, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 