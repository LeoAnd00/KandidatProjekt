# Creates plots for high and low nitrogen addition after starvation, with Sch9 as readout

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")

# Pre-shift => Glutamine_ext = 0, nitrogen starvation
include("../Model/ODE_methods.jl")
u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) 

# Low glutamine => Glutamine_ext = 0.3
tspan = (0.0, 30.0) # [min]

# Solve the ODEs under additon of low amounts of glutamine as nitrogen source
sol_low = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, p_var)

include("../../Data/exp_data_norm.jl")

# Creates the plot for low glutamine addition after starvation for Sch9
plot8 = scatter(
    t_Sch9P_glutamine_L, data_Sch9P_glutamine_L, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8
    )

plot!(sol_low, vars = Sch9, legend = false, color=color_jalihal, lw=2.0)
xlabel!("t [min]")
ylabel!("Sch9")
ylims!((0, 1.1))
xlims!((-0.5, last(tspan)*1.02))
title!("Kvävetillsättning (låg, 0.3)")
display(plot8)

file_name = "Sch9_nit_add_low_reconstr"
savefig(plot8, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot8, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 

# High Glutamine => Glutamine_ext = 1.0
# Solves the same ODEs with high amounts of glutamine 
sol_high = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var)

plot9 = scatter(
    t_Sch9P_glutamine_H, data_Sch9P_glutamine_H, markersize = 4.5, markercolor=RGB(0.35, 0.4, 1), markerstrokewidth=0.8
    )
plot!(sol_high, vars = Sch9, legend = false, color=color_jalihal, lw=2.0)
xlabel!("t [min]")
xlims!((-0.5, last(tspan)*1.02))
ylabel!("Sch9")
title!("Kvävetillsättning (hög, 1.0)")
display(plot9)

file_name = "Sch9_nit_add_high_reconstr"
savefig(plot9, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(plot9, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 