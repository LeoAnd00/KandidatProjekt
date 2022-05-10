# Creates plot for  Sch9 with a mutation gtr1_Delta under Glutamine addition

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")
include("../Model/ODE_methods.jl")

# Add the mutant gtr1_Delta by setting EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 
p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0

# Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1, nitogen starvation

# Returns u0 (initial conditions) from the pre-shift steady state
u0_SS = Steady_state_solver(p_const, p_var, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

# High glutamine => Glutamine_ext = 1, nitrogen addition
tspan = (0.0, 30.0) # [min]
# Solve the ODEs for the post-shift changes
sol = ODE_solver(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Plot the solution with experimental data from a separate file
plot7 = scatter(t_Sch9_gtr1Delta, data_Sch9_gtr1Delta, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8)

plot!(sol, vars=Sch9, legend=false, color=color_jalihal, lw=2.0)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
title!("Kvävetillsättning, gtr1\\Delta")
display(plot7)

file_name = "Sch9_nit_add_reconstr"
savefig(plot7, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot7, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 