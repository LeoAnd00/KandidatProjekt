# 

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0

# Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1 
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

# High glutamine => Glutamine_ext = 1
tspan = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")

plot1 = scatter(t_Sch9_gtr1Delta, data_Sch9_gtr1Delta)

plot!(sol, vars=Sch9, legend=false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
title!("Glutamine addition gtr1\\Delta, Glutamine_{ext} = 1.0")
display(plot1)

file_name = "Sch9_nit_add_reconstr"
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 