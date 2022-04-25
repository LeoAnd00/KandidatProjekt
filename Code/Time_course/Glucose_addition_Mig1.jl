# 

using DifferentialEquations
using ModelingToolkit
using Plots

default(titlefontsize=13)
default(dpi=300)

include("../Model/ODE_functions.jl") 

include("../Model/parameter_values.jl")

# Pre-shift => ATP, Carbon = 0
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

# Post-shift => ATP, Carbon = 1
tspan = (0.0, 20.0) # [min]
sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")
scatter(t_Mig1_glucose_relief, data_Mig1_glucose_relief, markersize = 5.5, markercolor="lightblue", markerstrokewidth=0.8)

plot1 = plot!(sol, vars=log(10, Mig1/(1e-5 + 1.0-Mig1)), color = "blue", lw=2.5, legend=false) # 1e-5 ensures denom != 0
xlabel!("t [min]")
ylabel!("log(Mig1/(Mig1_T-Mig1)")
ylims!((1.0, 1.6))
xticks!([0.0, 20.0])
# yticks!([1.2, 1.4])
title!("Glucose addition")
display(plot1)

file_name = "Mig1_gluc_add_reconstr"
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot1, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 
