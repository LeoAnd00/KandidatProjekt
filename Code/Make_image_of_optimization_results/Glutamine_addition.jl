# 

using DifferentialEquations
using ModelingToolkit
using Plots

function Glutamine_addition_L()

    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")
    include("../../Data/exp_data_norm.jl")
    u0_SS = Steady_state_solver(p_const, p_var_alt_1, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady 
    # state för parametrarna p

    # Low glutamine => Glutamine_ext = 0.3
    tspan = (0.0, 30.0) # [min]
    sol_low = ODE_solver(u0_SS, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, 
    p_var_alt_1)

    plot1 = scatter(t_Sch9P_glutamine_L, data_Sch9P_glutamine_L, 
    markersize = 5.5, color = color_data, markerstrokewidth=0.8, 
    labels = "Data")

    plot!(plot1, sol_low, vars = Sch9, color = Color_alt_1, lw=2.5, 
    labels= label_alt_1)


    u0_SS = Steady_state_solver(p_const, p_var_alt_2, (Carbon => 1.0, 
    ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady state för 
    # parametrarna p

    # Low glutamine => Glutamine_ext = 0.3
    tspan = (0.0, 30.0) # [min]
    sol_low = ODE_solver(u0_SS, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, 
    p_var_alt_2)

    plot!(plot1, sol_low, vars=Sch9, color = Color_alt_2, lw=2.5, 
    labels= label_alt_2)


    u0_SS = Steady_state_solver(p_const, p_var_jalihal, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady 
    # state för parametrarna p

    # Low glutamine => Glutamine_ext = 0.3
    tspan = (0.0, 30.0) # [min]
    sol_low = ODE_solver(u0_SS, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, 
    p_var_jalihal)

    plot!(plot1, sol_low, vars=Sch9, color = Color_jalihal, lw=2.5, 
    labels= label_jalihal)

    xlabel!("t [min]")
    ylabel!("Sch9")
    ylims!((0, 1.1))
    xlims!((-0.5, last(tspan)*1.02))
    title!("Glutamine addition, Low Glutamine (0.3)")
    return plot1
end

function Glutamine_addition_H()
    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")

    u0_SS = Steady_state_solver(p_const, p_var_alt_1, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady 
    # state för parametrarna p

    # Low glutamine => Glutamine_ext = 0.3
    tspan = (0.0, 30.0) # [min]
    sol_high = ODE_solver(u0_SS, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, 
    p_var_alt_1)

    plot2 = scatter(t_Sch9P_glutamine_H, data_Sch9P_glutamine_H, 
    markersize = 5.5, color = color_data, markerstrokewidth=0.8, 
    labels = "Data")

    plot!(plot2, sol_high, vars = Sch9, color = Color_alt_1, lw=2.5, 
    labels= label_alt_1)

    u0_SS = Steady_state_solver(p_const, p_var_alt_2, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady 
    # state för parametrarna p

    # Low glutamine => Glutamine_ext = 0.3
    tspan = (0.0, 30.0) # [min]
    sol_high = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, 
    Glutamine_ext => 1.0), tspan, p_const, p_var_alt_2)

    plot!(plot2, sol_high, vars=Sch9, color = Color_alt_2, lw=2.5, 
    labels= label_alt_2)

    u0_SS = Steady_state_solver(p_const, p_var_jalihal, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady 
    # state för parametrarna p

    # Low glutamine => Glutamine_ext = 0.3
    tspan = (0.0, 30.0) # [min]
    sol_high = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, 
    Glutamine_ext => 1.0), tspan, p_const, p_var_jalihal)

    plot!(plot2, sol_high, vars=Sch9, color = Color_jalihal, lw=2.5, 
    labels= label_jalihal)

    xlabel!("t [min]")
    xlims!((-0.5, last(tspan)*1.02))
    ylabel!("Sch9")
    title!("Glutamine addition, High glutamine")
    return plot2
end