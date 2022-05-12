# 

using DifferentialEquations
using ModelingToolkit
using Plots

function Glucose_starvation_Snf1()

    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")
    include("../../Data/exp_data_norm.jl")
    u0_SS = Steady_state_solver(p_const, p_var_alt_1, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady 
    # state för parametrarna p

    # Post-shift, Glucose starvation => Carbon, ATP = 0
    tspan_Snf1 = (0.0, 61.0) # [min]
    sol_Snf1 = ODE_solver(u0_SS, (Carbon => 0.0, ATP => 0.0, 
    Glutamine_ext => 1.0), tspan_Snf1, p_const, p_var_alt_1)

    plot1 = scatter(t_Snf1, data_Snf1, markersize = 5.5, 
    color = color_data, markerstrokewidth=0.8, label = "Data")

    plot!(plot1, sol_Snf1, vars=Snf1, color = Color_alt_1, lw=2.5, 
    labels= label_alt_1)


    u0_SS = Steady_state_solver(p_const, p_var_alt_2, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady 
    # state för parametrarna p

    # Post-shift, Glucose starvation => Carbon, ATP = 0
    tspan_Snf1 = (0.0, 61.0) # [min]
    sol_Snf1 = ODE_solver(u0_SS, 
    (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Snf1, p_const, 
    p_var_alt_2)

    plot!(plot1, sol_Snf1, vars=Snf1, color = Color_alt_2, lw=2.5, 
    labels= label_alt_2)

    u0_SS = Steady_state_solver(p_const, p_var_jalihal, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady 
    # state för parametrarna p

    # Post-shift, Glucose starvation => Carbon, ATP = 0
    tspan_Snf1 = (0.0, 61.0) # [min]
    sol_Snf1 = ODE_solver(u0_SS, 
    (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Snf1, p_const, 
    p_var_jalihal)

    plot!(plot1, sol_Snf1, vars=Snf1, color = Color_jalihal, lw=2.5, 
    labels= label_jalihal)

    xlabel!("t [min]")
    ylabel!("Snf1")
    ylims!((0.0, 1.1))
    xlims!((-0.5, last(tspan_Snf1)*1.02))
    title!("Glukossvältning Snf1")
    return plot1
end

function Glucose_starvation_Sch9()

    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")

    u0_SS = Steady_state_solver(p_const, p_var_alt_1, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) 
    tspan_Sch9 = (0.0, 30.0) # [min]
    sol_sch9 = ODE_solver(u0_SS, 
    (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, 
    p_var_alt_1)

    plot2 = scatter(t_Sch9_glucose_starve, data_Sch9_glucose_starve, 
    markersize = 5.5, color = color_data, markerstrokewidth=0.8, 
    label = "Data")

    plot!(plot2, sol_sch9, vars=Sch9, color = Color_alt_1, lw=2.5, 
    labels= label_alt_1)


    u0_SS = Steady_state_solver(p_const, p_var_alt_2, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) 
    tspan_Sch9 = (0.0, 30.0) # [min]
    sol_sch9 = ODE_solver(u0_SS, 
    (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, 
    p_var_alt_2)

    plot!(plot2, sol_sch9, vars=Sch9, color = Color_alt_2, lw=2.5, 
    labels= label_alt_2)

    u0_SS = Steady_state_solver(p_const, p_var_jalihal, 
    (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) 
    tspan_Sch9 = (0.0, 30.0) # [min]
    sol_sch9 = ODE_solver(u0_SS, 
    (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan_Sch9, p_const, 
    p_var_jalihal)

    plot!(plot2, sol_sch9, vars=Sch9, color = Color_jalihal, lw=2.5, 
    labels= label_jalihal)

    xlabel!("t [min]")
    ylabel!("Sch9")
    ylims!((0.0, 1.1))
    xlims!((-0.5, last(tspan_Sch9)*1.02))
    title!("Glukossvältning Sch9")
    return plot2
end
