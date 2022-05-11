# 

using DifferentialEquations
using ModelingToolkit
using Plots

function Glucose_addition_Sch9_cAMP()

    include("../Model/ODE_functions.jl") 

    include("../Model/parameter_values.jl")

    include("../Model/ODE_methods.jl")
    include("../../Data/exp_data_norm.jl")
    include("../../Results/Parameter_values/Results_from_optimization.jl")
    include("../../Data/exp_data_norm.jl")
    # Mutant Sch9_Delta => Sch9_T = 0
    p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0

    # Pre-shift => ATP, Carbon = 0
    u0_SS = Steady_state_solver(p_const, p_var_alt_1, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

    # Post-shift => ATP, Carbon = 1
    tspan = (0.0, 3.0) # [min]
    sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var_alt_1)

    plot1 = scatter(t_sch9Delta_cAMP, data_sch9Delta_cAMP, markersize = 5.5, color = color_data, markerstrokewidth=0.8, labels= "Data")

    plot!(plot1, sol, vars=cAMP, color = Color_alt_1, lw=2.5, labels= label_alt_1)


    u0_SS = Steady_state_solver(p_const, p_var_alt_2, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

    # Post-shift => ATP, Carbon = 1
    tspan = (0.0, 3.0) # [min]
    sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var_alt_2)

    plot!(plot1, sol, vars=cAMP, color = Color_alt_2, lw=2.5, labels= label_alt_2)

    u0_SS = Steady_state_solver(p_const, p_var_jalihal, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

    # Post-shift => ATP, Carbon = 1
    tspan = (0.0, 3.0) # [min]
    sol = ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var_jalihal)

    plot!(plot1, sol, vars=cAMP, color = Color_jalihal, lw=2.5, labels= label_jalihal)

    xlabel!("t [min]")
    ylabel!("cAMP")
    ylims!((0, 1.1))
    xlims!((-0.03, last(tspan)*1.01))
    title!("Glucose addition Sch9\\Delta")
    return plot1
end
