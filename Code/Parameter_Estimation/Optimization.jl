using LinearAlgebra
using DifferentialEquations
using Plots
using DiffEqSensitivity
using SteadyStateDiffEq
using ForwardDiff
using ProgressMeter
using ModelingToolkit
include("../Model/ODE_functions.jl") 
include("../Model/ODE_methods.jl")
include("../../Data/exp_data.jl")
include("Rearrange_the_gradient.jl")

"""
    main_while_loop(p_var_FWD, p_const_FWD, my, my2, abs_min_alpha, condition_num_to_use_SD)

    All calculations to be done per iteration.
"""
function main_while_loop(p_var_FWD, p_const_FWD, my, my2, abs_min_alpha, condition_num_to_use_SD)
    N = 0 
    
    #prog = ProgressUnknown(;dt=0.1,desc="Progress:", color=:green)
    step_dir = []
    step_length = []
    Hessian_inv = []
    grad = []
    prev_grad = []
    MK = []

    while true
        N += 1
        println("N:",N)
        #ProgressMeter.next!(prog)

        if N == 1
            try
                New_MK = calc_cost(p_var_FWD, p_const_FWD) # Calculates the cost of the initial parameter values
                MK = append!(MK, New_MK)
            catch
                return nothing
                break
            end
        end
        if N == 1
            grad = []
            try
                grad = calc_grad(p_var_FWD, p_const_FWD)
            catch
                break
            end
        end

        Hessian_inv = calc_hessian(prev_grad, grad, step_length, Hessian_inv, p_var_FWD)

        if cond(inv(Hessian_inv)) > condition_num_to_use_SD || N == 1 # Atm I'm not really using cond(inv(Hessian_inv)) > condition_num_to_use_SD
            println(" USING STEEPEST ")                                 # I'm doing steepest descent for the first iteration, N == 1
            step_dir = -grad
            step_length = steepest_descent_calc_length(step_dir, p_var_FWD, p_const_FWD, abs_min_alpha, MK)  
        else
            step_dir = calc_direction(grad, Hessian_inv)
            step_length, alpha_control = calc_length(step_dir, p_var_FWD, p_const_FWD, MK, grad, my, my2, abs_min_alpha)
            if alpha_control == 0
                println(" USING STEEPEST ")
                step_dir = -grad
                step_length = steepest_descent_calc_length(step_dir, p_var_FWD, p_const_FWD, abs_min_alpha, MK) 
            end
        end

        p_values_process = exp.(log.(last.(p_var_FWD))+step_length) # new p values based on step length calculated previously

        p_values_dual_temp = [Pair(first.(p_var_FWD)[1], p_values_process'[1])]
        for i in range(2, length(last.(p_var_FWD)))
            B = Pair(first.(p_var_FWD)[i], p_values_process'[i])
            p_values_dual_temp = vcat(p_values_dual_temp,B)
        end
        p_var_FWD = p_values_dual_temp # Makes it into a Pair vector

        New_MK = calc_cost(p_var_FWD, p_const_FWD)

        MK = append!(MK, New_MK)

        prev_grad = grad

        grad = []
        try
            grad = calc_grad(p_var_FWD, p_const_FWD)
        catch
            break
        end

        Number_of_fulfilled_criteria = termination_criteria(MK, grad, step_length, p_var_FWD)
        if Number_of_fulfilled_criteria >= 2
            break
        end
    end
    #ProgressMeter.finish!(prog)

    return MK, p_var_FWD, N
end

"""
    termination_criteria(MK, grad, step_length, p_var_FWD)

    Checks how many termination criterias are fulfilled.
"""
function termination_criteria(MK, grad, step_length, p_var_FWD)

    epsilon_1 = 1e-6
    epsilon_2 = 1e-6
    epsilon_3 = 1e-6

    Number_of_fulfilled_criteria = 0
    #criteria_1
    if norm(grad) <= epsilon_1*(1+abs(MK[end]))
        Number_of_fulfilled_criteria += 1
    end
    #criteria_2
    if MK[end-1] - MK[end] <= epsilon_2*(1+abs(MK[end]))
        Number_of_fulfilled_criteria += 1
    end
    #criteria_3
    if norm(step_length) <= epsilon_3*(1+norm(last.(p_var_FWD)))
        Number_of_fulfilled_criteria += 1
    end

    return Number_of_fulfilled_criteria
end

"""
    calc_hessian(prev_grad, grad, step_length, Hessian_inv, p_var_FWD)

    Calculate the inverse of the hessian.
"""
function calc_hessian(prev_grad, grad, step_length, Hessian_inv, p_var_FWD)

    if prev_grad == []
        dimentions = size(last.(p_var_FWD))[1]
        Hessian = Matrix(I,dimentions,dimentions)
        Hessian_inv  = inv(Hessian)
    else
        y = grad .- prev_grad
        #Hessian invers with "Broyden–Fletcher–Goldfarb–Shanno algorithm"
        Hessian_inv = Hessian_inv .+ (step_length'*y + y'*(Hessian_inv*y))*step_length*step_length'/((y'*step_length)^2) .- (Hessian_inv*y*step_length' .+ step_length*y'*Hessian_inv)/(step_length'*y)
    end
    
    return Hessian_inv
end

"""
    calc_direction(grad, Hessian_inv)

    Calculate the step direction.
"""
function calc_direction(grad, Hessian_inv)

    step_dir = -Hessian_inv*grad

    return step_dir
end

"""
    calc_length(step_dir, p_var_FWD, p_const_FWD, MK, grad, my, my2, abs_min_alpha)

    Finds a step length that fulfills the Armijo and Wolfe condition
"""
function calc_length(step_dir, p_var_FWD, p_const_FWD, MK, grad, my, my2, abs_min_alpha)
    #gå i stegriktningen
    alpha = 1
    p_var_new = []

    while true
        #Minimum alpha allowed so it doesn't go to inf.
        if alpha < abs_min_alpha 
            println(" BAD ALPHA ")
            alpha = 0 # used to switch to steepest descent if alpha < abs_min_alpha for quasi-newton
            break
        end

        p_values_prev = exp.(log.(last.(p_var_FWD))+alpha*step_dir)

        p_values_dual_temp = [Pair(first.(p_var_FWD)[1], p_values_prev'[1])]
        for i in range(2, length(last.(p_var_FWD)))
            B = Pair(first.(p_var_FWD)[i], p_values_prev'[i])
            p_values_dual_temp = vcat(p_values_dual_temp,B)
        end
        p_var_new = p_values_dual_temp

        MK_comp_val = 0
        try
            MK_comp_val = calc_cost(p_var_new, p_const_FWD) # Checks if the cost can be calculated for the parameters
        catch 
            println(" Can't calculate MK in Quasi, alpha = ",alpha)
            alpha  = alpha/2
            continue
        end
        println(" QUASI_MK: ", MK_comp_val)

        MK_diff = MK_comp_val - MK[end]

        Armijo_cond = my * alpha * step_dir[:,end]' * grad[:,end]


        if MK_diff > Armijo_cond
            alpha  = alpha/2
        else
            println("Fullfilled first condition!")
            D = 0
            try
                D = calc_grad(p_var_new, p_const_FWD)
            catch 
                println(" Can't calculate grad in Quasi, alpha = ",alpha)
                alpha  = alpha/2
                continue
            end

            step_direct_grad = -step_dir[:,end]' * D

            Wolfe_cond = -my2 * step_dir[:,end]' * grad[:,end]

            if abs(step_direct_grad) < abs(Wolfe_cond) #strong Wolfe conditions
                break
            else
                alpha = alpha/2
            end
        end

    end

    step_length = alpha*step_dir[:,end]

    return step_length, alpha
end

"""
    calc_grad(p_var_FWD, p_const_FWD)

    Calculate the gradient of the cost function.
"""
function calc_grad(p_var_FWD, p_const_FWD)

    function dMK_dp_Glucose_addition_Mig1(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)

            u_approx = getindex.(sol.u, 20) #ODE for MIG1 on row 20 in ODE system

            MK_for_diff = sum((data_for_ode - (log.(10,u_approx./(1e-5 .+ 1.0 .- u_approx)))).^2)
            return MK_for_diff
        end
    
        data_for_ode = raw_data_Mig1_glucose_relief
        timelist_for_ode = t_Mig1_glucose_relief
        tspan = (0.0, 20.0) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)
        #u0_SS = Steady_state_solver(p_var_FWD, p_const_FWD, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) 

        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Glucose_addition_Sch9(p_var_FWD,p_const_FWD)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 6)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_sch9Delta_cAMP, 0.052, 1.227)
        timelist_for_ode = t_sch9Delta_cAMP
        tspan = (0.0, 3.5) # [min]
        p_const_FWD[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0
        p_values = vcat(p_var_FWD,p_const_FWD)
        u0_SS = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))
        p_const_FWD[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 1.0

        return dMK_dp
    end

    function dMK_dp_Glucose_addition_cAMP(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 6) 

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_cAMP, 0.052, 1.227)
        timelist_for_ode = t_cAMP
        tspan = (0.0, 3.5) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Glucose_addition_Sch9_p(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12) 

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Sch9_glucose_relief)
        timelist_for_ode = t_Sch9_glucose_relief
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Glucose_starvation_Snf1_p(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 10)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Snf1)
        timelist_for_ode = t_Snf1
        tspan = (0.0, 62) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Glucose_starvation_Sch9_p(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Sch9_glucose_starve)
        timelist_for_ode = t_Sch9_glucose_starve
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Sch9P_glutamine_L(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Sch9P_glutamine_L, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_L
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Sch9P_glutamine_H(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Sch9P_glutamine_H, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_H
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var_FWD,p_const_FWD)

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))

        return dMK_dp
    end

    function dMK_dp_Glutamine_addition_Sch9_gtr1Delta(p_var_FWD,p_const_FWD)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Sch9_gtr1Delta, 4.82, 52.27)
        timelist_for_ode = t_Sch9_gtr1Delta
        tspan = (0.0, 32) # [min]
        w_torc_ego_true = p_var_FWD[Get_index(p_var_lookup_table, "w_torc_ego")]
        w_torc_egoin_true = p_var_FWD[Get_index(p_var_lookup_table, "w_torc_egoin")]

        # Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
        p_const_FWD[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0
        p_values = vcat(p_var_FWD,p_const_FWD)
        u0_SS = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 

        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))
        p_const_FWD[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 1.0
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego_true
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin_true

        return dMK_dp
    end

    function dMK_dp_Rapamycin_treatment(p_var_FWD,p_const_FWD, u0_SS)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 25)

            MK_for_diff = sum((data_for_ode - (u_approx./(1e-3 .+ u0_SS[25]))).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(raw_data_Rib_rap)
        timelist_for_ode = t_Rib_rap
        tspan = (0.0, 92.0) # [min]

        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
        p_const_FWD[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0
        p_values = vcat(p_var_FWD,p_const_FWD)
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))
        p_const_FWD[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 1.0

        return dMK_dp
    end

    u0_SS_A_C_G = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
    #println("u0_SS_A_C_G:",u0_SS_A_C_G)
    dMK_dp_Glucose_starvation_Snf1_p_solution = dMK_dp_Glucose_starvation_Snf1_p(p_var_FWD,p_const_FWD, u0_SS_A_C_G)
    dMK_dp_Glucose_starvation_Sch9_p_solution = dMK_dp_Glucose_starvation_Sch9_p(p_var_FWD,p_const_FWD, u0_SS_A_C_G)
    dMK_dp_Rapamycin_treatment_solution = dMK_dp_Rapamycin_treatment(p_var_FWD,p_const_FWD, u0_SS_A_C_G)

    u0_SS_G = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0))
    #println("u0_SS_G:",u0_SS_G)
    dMK_dp_Glucose_addition_Mig1_solution = dMK_dp_Glucose_addition_Mig1(p_var_FWD,p_const_FWD, u0_SS_G)
    dMK_dp_Glucose_addition_cAMP_solution = dMK_dp_Glucose_addition_cAMP(p_var_FWD,p_const_FWD, u0_SS_G)
    dMK_dp_Glucose_addition_Sch9_p_solution = dMK_dp_Glucose_addition_Sch9_p(p_var_FWD,p_const_FWD, u0_SS_G)

    u0_SS_A_C = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0))
    #println("u0_SS_A_C:",u0_SS_A_C)
    dMK_dp_Sch9P_glutamine_L_solution = dMK_dp_Sch9P_glutamine_L(p_var_FWD,p_const_FWD, u0_SS_A_C)
    dMK_dp_Sch9P_glutamine_H_solution = dMK_dp_Sch9P_glutamine_H(p_var_FWD,p_const_FWD, u0_SS_A_C)

    dMK_dp_Glucose_addition_Sch9_solution = dMK_dp_Glucose_addition_Sch9(p_var_FWD,p_const_FWD)
    dMK_dp_Glutamine_addition_Sch9_gtr1Delta_solution = dMK_dp_Glutamine_addition_Sch9_gtr1Delta(p_var_FWD,p_const_FWD)

    True_grad = dMK_dp_Glucose_addition_Mig1_solution + dMK_dp_Glucose_addition_Sch9_solution + dMK_dp_Glucose_addition_cAMP_solution + dMK_dp_Glucose_addition_Sch9_p_solution + dMK_dp_Glucose_starvation_Snf1_p_solution + dMK_dp_Glucose_starvation_Sch9_p_solution + dMK_dp_Sch9P_glutamine_L_solution + dMK_dp_Sch9P_glutamine_H_solution + dMK_dp_Glutamine_addition_Sch9_gtr1Delta_solution + dMK_dp_Rapamycin_treatment_solution

    True_grad = Rearrange_the_gradient(p_var_FWD, p_const_FWD, True_grad)

    return True_grad
end

"""
    calc_cost(p_var_FWD, p_const_FWD)

    Calculate the sum of least square / cost for all chosen ODE:s
"""
function calc_cost(p_var_FWD, p_const_FWD)

    function Calc_cost_Glucose_addition_Mig1(p_var_FWD, p_const_FWD, u0_SS)

        data_for_ode = raw_data_Mig1_glucose_relief
        timelist_for_ode = t_Mig1_glucose_relief
        tspan = (0.0, 20.0) # [min]

        #u0_SS = Steady_state_solver(p_const, p_var, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

        sol = ODE_solver_FWD(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
        
        u_approx = getindex.(sol.u, 20) #ODE for MIG1 on row 20 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glucose_addition_Mig1 = sum((data_for_ode - (log.(10,u_approx./(1e-5 .+ 1.0 .- u_approx)))).^2)
    
        return MK_Glucose_addition_Mig1
    end

    function Calc_cost_Glucose_addition_Sch9(p_var_FWD, p_const_FWD)

        data_for_ode = Minmaxnorm(raw_data_sch9Delta_cAMP, 0.052, 1.227)
        timelist_for_ode = t_sch9Delta_cAMP
        tspan = (0.0, 3.5) # [min]
        # Mutant Sch9_Delta => Sch9_T = 0
        p_const_FWD[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0
    
        # Pre-shift => ATP, Carbon = 0
        u0_SS = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift => ATP, Carbon = 1
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
        
        u_approx = getindex.(sol.u, 6) #ODE for cAMP on row 6 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glucose_addition_Sch9 = sum((data_for_ode - u_approx).^2)

        p_const_FWD[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 1.0
    
        return MK_Glucose_addition_Sch9
    end

    function Calc_cost_Glucose_addition_cAMP(p_var_FWD, p_const_FWD, u0_SS)

        data_for_ode = Minmaxnorm(raw_data_cAMP, 0.052, 1.227)
        timelist_for_ode = t_cAMP
        tspan = (0.0, 3.5) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p

        # Post-shift, Glucose addition => ATP, Carbon = 1
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 6) #ODE for cAMP on row 6 in ODE system
            
        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glucose_addition_cAMP = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_addition_cAMP
    end

    function Calc_cost_Glucose_addition_Sch9_p(p_var_FWD, p_const_FWD, u0_SS)

        data_for_ode = Minmaxnorm(raw_data_Sch9_glucose_relief)
        timelist_for_ode = t_Sch9_glucose_relief
        tspan = (0.0, 30) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glucose_addition_Sch9_p = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_addition_Sch9_p
    end

    function Calc_cost_Glucose_starvation_Snf1_p(p_var_FWD, p_const_FWD,u0_SS)

        data_for_ode = Minmaxnorm(raw_data_Snf1)
        timelist_for_ode = t_Snf1
        tspan = (0.0, 62) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift, Glucose starvation => Carbon, ATP = 0
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
        u_approx = getindex.(sol.u, 10) #ODE for Snf1 on row 10 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glucose_starvation_Snf1_p = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_starvation_Snf1_p
    end

    function Calc_cost_Glucose_starvation_Sch9_p(p_var_FWD, p_const_FWD,u0_SS)

        data_for_ode = Minmaxnorm(raw_data_Sch9_glucose_starve)
        timelist_for_ode = t_Sch9_glucose_starve
        tspan = (0.0, 30) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glucose_starvation_Sch9_p = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_starvation_Sch9_p
    end

    function Calc_cost_Sch9P_glutamine_L(p_var_FWD, p_const_FWD, u0_SS)

        data_for_ode = Minmaxnorm(raw_data_Sch9P_glutamine_L, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_L
        tspan = (0.0, 32) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady state för parametrarna p
    
        # Low glutamine => Glutamine_ext = 0.3
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Sch9P_glutamine_L = sum((data_for_ode - u_approx).^2)

        return MK_Sch9P_glutamine_L
    end

    function Calc_cost_Sch9P_glutamine_H(p_var_FWD, p_const_FWD, u0_SS)

        data_for_ode = Minmaxnorm(raw_data_Sch9P_glutamine_H, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_H
        tspan = (0.0, 32) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady state för parametrarna p
    
        # High Glutamine => Glutamine_ext = 1.0
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system
 
        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Sch9P_glutamine_H = sum((data_for_ode - u_approx).^2)

        return MK_Sch9P_glutamine_H
    end

    function Calc_cost_Glutamine_addition_Sch9_gtr1Delta(p_var_FWD, p_const_FWD)

        data_for_ode = Minmaxnorm(raw_data_Sch9_gtr1Delta, 4.82, 52.27)
        timelist_for_ode = t_Sch9_gtr1Delta
        tspan = (0.0, 32) # [min]
    
        w_torc_ego_true = p_var_FWD[Get_index(p_var_lookup_table, "w_torc_ego")]
        w_torc_egoin_true = p_var_FWD[Get_index(p_var_lookup_table, "w_torc_egoin")]
    
        # Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
        p_const_FWD[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0
    
        # Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1 
        u0_SS = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 
    
        # High glutamine => Glutamine_ext = 1
        sol = ODE_solver_FWD(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system
           
        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Glutamine_addition_Sch9_gtr1Delta = sum((data_for_ode - u_approx).^2)
    
        p_const_FWD[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 1.0
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego_true
        p_var_FWD[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin_true
    
        return MK_Glutamine_addition_Sch9_gtr1Delta
    end

    function Calc_cost_Rapamycin_treatment(p_var_FWD, p_const_FWD, u0_SS)

        data_for_ode = Minmaxnorm(raw_data_Rib_rap)
        timelist_for_ode = t_Rib_rap
        tspan = (0.0, 92.0) # [min]
    
        #u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift: Rapamycin treatment => TORC1_T = 0.0 
        p_const_FWD[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const_FWD, p_var_FWD, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 25) #ODE for RIB on row 25 in ODE system

        #Check if sol was a success
        u_approx = Check_if_cost_is_a_success(sol, u_approx)

        MK_Rapamycin_treatment = sum((data_for_ode - (u_approx./(1e-3 + u0_SS[25]))).^2)
    
        p_const_FWD[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 1.0
    
        return MK_Rapamycin_treatment
    end
    
    u0_SS_A_C_G = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
    MK_Glucose_starvation_Snf1_p = Calc_cost_Glucose_starvation_Snf1_p(p_var_FWD, p_const_FWD, u0_SS_A_C_G)
    MK_Glucose_starvation_Sch9_p = Calc_cost_Glucose_starvation_Sch9_p(p_var_FWD, p_const_FWD, u0_SS_A_C_G)
    MK_Rapamycin_treatment = Calc_cost_Rapamycin_treatment(p_var_FWD, p_const_FWD, u0_SS_A_C_G)

    u0_SS_G = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0))
    MK_Glucose_addition_Mig1 = Calc_cost_Glucose_addition_Mig1(p_var_FWD, p_const_FWD, u0_SS_G)
    MK_Glucose_addition_cAMP = Calc_cost_Glucose_addition_cAMP(p_var_FWD, p_const_FWD, u0_SS_G)
    MK_Glucose_addition_Sch9_p = Calc_cost_Glucose_addition_Sch9_p(p_var_FWD, p_const_FWD, u0_SS_G)

    u0_SS_A_C = Steady_state_solver_FWD(p_const_FWD, p_var_FWD, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0))
    MK_Sch9P_glutamine_L = Calc_cost_Sch9P_glutamine_L(p_var_FWD, p_const_FWD, u0_SS_A_C)
    MK_Sch9P_glutamine_H = Calc_cost_Sch9P_glutamine_H(p_var_FWD, p_const_FWD, u0_SS_A_C)

    MK_Glucose_addition_Sch9 = Calc_cost_Glucose_addition_Sch9(p_var_FWD, p_const_FWD)
    MK_Glutamine_addition_Sch9_gtr1Delta = Calc_cost_Glutamine_addition_Sch9_gtr1Delta(p_var_FWD, p_const_FWD)

    MK = MK_Glucose_addition_Mig1 + MK_Glucose_addition_Sch9 + MK_Glucose_addition_cAMP + MK_Glucose_addition_Sch9_p + MK_Glucose_starvation_Snf1_p + MK_Glucose_starvation_Sch9_p + MK_Sch9P_glutamine_L + MK_Sch9P_glutamine_H + MK_Glutamine_addition_Sch9_gtr1Delta + MK_Rapamycin_treatment

    return MK
end

"""
    Check_if_cost_is_a_success(sol, u_approx)

    Checks if sol.retcode != :Success
"""
function Check_if_cost_is_a_success(sol, u_approx)
    if sol.retcode != :Success
        println(" Can't Calc Cost ")
        u_approx = "Not working"
    end
    return u_approx
end

"""
    plot_results(Bild, Bild_format, MK, Number_of_iterations)

    Plot the chosen ODE with or without data points.
"""
function plot_results(Bild, Bild_format, MK, Number_of_iterations)
    if Bild == true
        gr()
        plotd = plot((0:Number_of_iterations),MK,legend=:topright)
        plot!(plotd,(0:Number_of_iterations),MK, seriestype = :scatter)
        xlabel!("Iteration")
        ylabel!("Kostnad")
        title!("Kostnad Vs Iteration")
        if Bild_format == "PNG"
            savefig(plotd, pwd()*"/Results/MK_plot.png")
        elseif Bild_format == "SVG"
            savefig(plotd, pwd()*"/Results/MK_plot.svg")
        end
        display(plotd)
    end
end

"""
    steepest_descent_calc_length(step_dir, p_var_FWD, p_const_FWD, abs_min_alpha, MK)  

    Finds a step length for steepest descent.
"""
function steepest_descent_calc_length(step_dir, p_var_FWD, p_const_FWD, abs_min_alpha, MK)  
    alpha0 = 1e-9
    alpha = 1e-9
    p_var_temp = []
    abs_min_alpha = abs_min_alpha*1e-10

    while true

        p_values_prev = exp.(log.(last.(p_var_FWD)) + (alpha-alpha0)*step_dir)

        println(" Alpha: ", alpha)

        p_values_dual_temp = [Pair(first.(p_var_FWD)[1], p_values_prev'[1])]
        for i in range(2, length(last.(p_var_FWD)))
            B = Pair(first.(p_var_FWD)[i], p_values_prev'[i])
            p_values_dual_temp = vcat(p_values_dual_temp,B)
        end
        p_var_temp = p_values_dual_temp

        prev_MK = 0
        try
            prev_MK = calc_cost(p_var_temp, p_const_FWD)
        catch 
            if alpha <= abs_min_alpha
                alpha = 0.0
                break
            end
            if alpha > alpha0
                alpha = alpha-alpha0
                break
            else
                alpha0 = alpha0/10
                alpha = alpha0
                continue
            end
        end

        p_values_new = exp.(log.(last.(p_var_FWD))+alpha*step_dir)

        p_values_dual_temp = [Pair(first.(p_var_FWD)[1], p_values_new'[1])]
        for i in range(2, length(last.(p_var_FWD)))
            B = Pair(first.(p_var_FWD)[i], p_values_new'[i])
            p_values_dual_temp = vcat(p_values_dual_temp,B)
        end

        p_var_new = p_values_dual_temp

        new_MK = 0
        try
            new_MK = calc_cost(p_var_new, p_const_FWD)
        catch
            if alpha <= abs_min_alpha
                alpha = 0.0
                break
            end
            if alpha > alpha0
                alpha = alpha-alpha0
                break
            else
                alpha0 = alpha0/10
                alpha = alpha0
                continue
            end
        end
        println("New_MK: ",new_MK)


        # If alpha = 0 it will try reducing alpha0 until alpha is not 0 or alpha0 > abs_min_alpha
        if new_MK > prev_MK
            alpha = alpha-alpha0
            if alpha == 0
                if alpha0 <= abs_min_alpha
                    break
                end
                alpha0 = alpha0/10
                alpha = alpha0
                continue
            end
            break
        else
            if alpha <= abs_min_alpha
                alpha = 0.0
                break
            end
            if alpha >= 10*alpha0
                alpha0 = 10*alpha0
            end
            
            if round(alpha) >= 1e6
                break
            end

            alpha = alpha + alpha0
        end
    end

    step_length =alpha*step_dir[:,end]

    println(" Final_ALPHA: ",alpha)

    return step_length
end

"""
    main_FWD_optimization(p_0_var)

    Here you choose different options and values for the entire optimisation process.
"""
function main_FWD_optimization(p_0_var)
    println("Works")
    include("../Model/parameter_values.jl")
    println("works2")
    p_var = p_0_var
    condition_num_to_use_SD = 1e100 #Not using atm, hence a big number 1e100

    #calc_length options:
    my = 1e-4
    my2 = 0.95 
    #my2 = 1e3
    abs_min_alpha = 1e-8

    ####################Quasi-Newton######################
    MK, p_var_reslut, Number_of_iterations = @time main_while_loop(p_var, p_const, my, my2, abs_min_alpha, condition_num_to_use_SD)
    #######################################################

    #plot cost vs iterations
    Bild = false
    Bild_format = "PNG" #Defined for PNG and SVG
    plot_results(Bild, Bild_format, MK, Number_of_iterations)
    return MK, p_var_reslut
end

#include("../Model/parameter_values.jl")
#p_0_var = p_var
#MK, p_var_result = main_FWD_optimization(p_0_var)
#println(p_var_result)
#println(MK)
