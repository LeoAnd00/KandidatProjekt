# 

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift, Glucose addition => ATP, Carbon = 1
tspan_cAMP = (0.0, 3.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_cAMP, p_const, p_var)

include("../../Data/exp_data_norm.jl")

plot1 = scatter(t_cAMP, data_cAMP)

plot!(sol, vars = cAMP, legend = false)
xlabel!("t [min]")
ylabel!("cAMP")
ylims!((0.0, 1.02))
xlims!((-0.02, last(tspan_cAMP)*1.02))
title!("Glucose addition")
display(plot1)

file_name = "cAMP_gluc_add_reconstr"
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 

tspan_Sch9 = (0.0, 30.0) # [min]
sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, p_var)

plot2 = scatter(t_Sch9_glucose_relief, data_Sch9_glucose_relief)

plot!(sol, vars = Sch9, legend = false)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0.0, 1.02))
xlims!((-0.3, last(tspan_Sch9)*1.02))
title!("Glucose addition")
display(plot2)

file_name = "Snf1_gluc_starve_reconstr"
savefig(plot2, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot2, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 
