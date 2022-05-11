# 

using DifferentialEquations
using ModelingToolkit
using Plots

function Glutamine_addition_gtr1()

    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")
    include("../../Data/exp_data_norm.jl")

    # Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
    p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
    p_var_alt_1[Get_index(p_var_lookup_table, "w_torc_ego")] = 
    w_torc_ego => 0.0
    p_var_alt_1[Get_index(p_var_lookup_table, "w_torc_egoin")] = 
    w_torc_egoin => 0.0

    # Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1 

    u0_SS = Steady_state_solver(p_const, p_var_alt_1, 
    (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

    # High glutamine => Glutamine_ext = 1
    tspan = (0.0, 30.0) # [min]
    sol = ODE_solver(u0_SS, 
    (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const, 
    p_var_alt_1)

    plot1 = scatter(t_Sch9_gtr1Delta, data_Sch9_gtr1Delta, markersize = 5.5, 
    color = color_data, markerstrokewidth=0.8, label = "Data")

    plot!(plot1, sol, vars=Sch9, color = Color_alt_1, lw=2.5, 
    labels= label_alt_1)

    p_var_alt_2[Get_index(p_var_lookup_table, "w_torc_ego")] = 
    w_torc_ego => 0.0
    p_var_alt_2[Get_index(p_var_lookup_table, "w_torc_egoin")] = 
    w_torc_egoin => 0.0
    u0_SS = Steady_state_solver(p_const, p_var_alt_2, 
    (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

    # High glutamine => Glutamine_ext = 1
    tspan = (0.0, 30.0) # [min]
    sol = ODE_solver(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0),
    tspan, p_const, p_var_alt_2)

    plot!(plot1, sol, vars=Sch9, color = Color_alt_2, lw=2.5, 
    labels= label_alt_2)


    p_var_jalihal[Get_index(p_var_lookup_table, "w_torc_ego")] = 
    w_torc_ego => 0.0
    p_var_jalihal[Get_index(p_var_lookup_table, "w_torc_egoin")] = 
    w_torc_egoin => 0.0
    u0_SS = Steady_state_solver(p_const, p_var_jalihal, 
    (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

    # High glutamine => Glutamine_ext = 1
    tspan = (0.0, 30.0) # [min]
    sol = ODE_solver(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0),
    tspan, p_const, p_var_jalihal)

    plot!(plot1, sol, vars=Sch9, color = Color_jalihal, lw=2.5, 
    labels= label_jalihal)


    xlabel!("t [min]")
    ylabel!("Sch9")
    ylims!((0.0, 1.1))
    xlims!((-0.5, last(tspan)*1.02))
    title!("Glutamine addition gtr1\\Delta, Glutamine_{ext} = 1.0")
    return plot1
end