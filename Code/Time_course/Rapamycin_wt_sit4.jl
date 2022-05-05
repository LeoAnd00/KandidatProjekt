# Creates the plot for adding Rapamycin treatment after a specifice time  

using DifferentialEquations
using ModelingToolkit
using Plots

default(dpi = 300)
default(titlefontsize=13)

include("../Model/ODE_functions.jl") 

########## Sit4- mutant #########
include("../Model/parameter_values.jl")

# Pre-shift Sit4 mutant which means that TORC1_T = 0
p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0

include("../Model/ODE_methods.jl")
u0_SS_sit4 = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returns steady state parameters p

# Post-shift, Rapamycin treatment Sit4 => TORC1_T = 0 => no change in parameter values
tspan_sit4 = (0.0, 150.0) # [min]
sol_sit4 = ODE_solver(u0_SS_sit4, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_sit4, p_const, p_var) # SOlves the ODEs

########## wild type #########
include("../Model/parameter_values.jl")

# Pre-shift wild typ => no change in parameter values

u0_SS_wt = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returns steady state for parameters p

# Post-shift, Rapamycin treatment wild type => TORC1_T = 0
p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0

tspan_wt = (60.0, 150.0)
sol_wt = ODE_solver(u0_SS_wt, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan_wt, p_const, p_var) # Solves ODEs

# Plot the results 
plot(sol_sit4, vars = 1-Gln3, label = "Sit4\\Delta", lw = 4, ls=:dash, legend=:right)
plot!(sol_wt, vars = Gln3, label = "vildtyp", lw = 4, ls=:dash, color = "red")
p1 = plot!([0.0, 60.0], [u0_SS_wt[16], u0_SS_wt[16]], ls=:dash, lw=4, color = "red", label="")
vline!([60], color = "gray", label = "")
# annotate!(65, 0.95, "Rapamycin treatment")

xlabel!("t [min]")
ylabel!("Gln3") 
xlims!((0, 150))
ylims!((0, 1.0))
title!("Rapamycinbehandling")
display(p1)

file_name = "Gln3_rap_treat_wt_sit4_reconstr"
savefig(p1, pwd()*"/Results/Time_course_reconstructed/pixel_images/"*file_name*".png")    
savefig(p1, pwd()*"/Results/Time_course_reconstructed/vector_images/"*file_name*".pdf") 