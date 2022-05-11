# 

using DifferentialEquations
using ModelingToolkit
using Plots

function Glucose_addition_Mig1()
    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    # Pre-shift => ATP, Carbon = 0
    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")
    include("../../Data/exp_data_norm.jl")

    plot1 = scatter(t_Mig1_glucose_relief, data_Mig1_glucose_relief, 
    markersize = 5.5, color = color_data, markerstrokewidth=0.8, 
    labels= "Data")

    u0_SS = Steady_state_solver(p_const, p_var_alt_1, (ATP => 0.0, 
    Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för 
    # parametrarna p

    # Post-shift => ATP, Carbon = 1
    tspan = (0.0, 20.0) # [min]
    sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0),
    tspan, p_const, p_var_alt_1)

    plot!(plot1, sol, vars=log(10, Mig1/(1e-5 + 1.0-Mig1)), 
    color = Color_alt_1, lw=2.5, labels= label_alt_1)

    u0_SS = Steady_state_solver(p_const, p_var_alt_2, (ATP => 0.0, 
    Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för 
    # parametrarna p

    # Post-shift => ATP, Carbon = 1
    tspan = (0.0, 20.0) # [min]
    sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0),
    tspan, p_const, p_var_alt_2)

    plot!(plot1, sol, vars=log(10, Mig1/(1e-5 + 1.0-Mig1)), 
    color = Color_alt_2, lw=2.5, labels= label_alt_2)


    u0_SS = Steady_state_solver(p_const, p_var_jalihal, 
    (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady 
    # state för parametrarna p

    # Post-shift => ATP, Carbon = 1
    tspan = (0.0, 20.0) # [min]
    sol = ODE_solver(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0),
    tspan, p_const, p_var_jalihal)

    plot!(plot1,sol, vars=log(10, Mig1/(1e-5 + 1.0-Mig1)), 
    color = Color_jalihal, lw=2.5, labels= label_jalihal) 

    xlabel!("t [min]")
    ylabel!("log(Mig1/(Mig1_T-Mig1)")
    ylims!((0.9, 1.6))
    title!("Glucose addition")

    return plot1
end


