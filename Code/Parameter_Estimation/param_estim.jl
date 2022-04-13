using DifferentialEquations
using LinearAlgebra
using ForwardDiff

"""
    est_p_var(ODE_system, tvals, initial_values, start_guess[, n_parameters])

    Estimate parameters for ODE_system.

    If n_parameters is unspecified, estimation is done for all parameters in start_guess.
"""
function est_p_var(ODE_sys, tvals, init_vals, p_const, p_var, pre_shift, post_shift)

    tspan = [first(tvals),last(tvals)]

    #Information updated in each iteration.
    prev_p_vars = [] #list of all previous parameter values, "memory"
    prev_costs = [] #list of corresponding costs, "memory"
    prev_H = [1 0; 0 1] #last Hessian approximation
    prev_grad = [0,0] #last gradient

    #Constants: epsilon used in termination criteria and mu, nu used in step length calculations.
    epsilon,mu,nu = initiate_constants()

    #---------------------#---------------------
    #Cost as function of parameters. Uses cost_fun(pred) below.
    function cost_fun_of_p(p_var)
        return calc_cost(p_var, p_const)
    end

    #---------------------#---------------------
    #Conditions used when calculating the step length alpha.

    #General condition on alpha: positive parameters.
    function GC(current_p_vars,search_dir,alpha)
        for i in 1:n_ps
            if current_p_vars[i]+alpha*search_dir[i] < 0
                return false
            end
        end
        return true
    end

    #Condition on alpha used for steepest descent: new cost lower than previous cost.
    function SDC(current_p_vars,search_dir,alpha)
        return cost_fun_of_p(current_p_vars+alpha*search_dir) <= cost_fun_of_p(current_p_vars)
    end

    #Wolfe condition 1 (aka Armijo condition).
    function WC1(current_p_vars,grad,search_dir,alpha)
        return cost_fun_of_p(current_p_vars+alpha*search_dir)-cost_fun_of_p(current_p_vars) <= mu*alpha*search_dir'*grad
    end

    #Wolfe condition 2 (strong form).
    function WC2(current_p_vars,grad,search_dir,alpha)
        return abs(search_dir'*calc_grad(exp_ODE_sys,init_vals,tspan,current_p_vars+alpha*search_dir,tvals)) >= abs(nu*search_dir'*grad)
    end

    #---------------------#---------------------
    #Functions to calculate step length (alpha).

    #For steepest descent.
    function steepest_descent_step(current_p_vars,search_dir)
        alpha = 1
        while !(SDC(current_p_vars,search_dir,alpha) && GC(current_p_vars,search_dir,alpha))
            alpha *= 1/2
        end
        return alpha
    end

    #Armijo step length rule (using Wolfe conditions).
    function armijo(current_p_vars,grad,search_dir)
        alpha = 1
        accept = true
        while !(WC1(current_p_vars,grad,search_dir,alpha) && WC2(current_p_vars,grad,search_dir,alpha) && GC(current_p_vars,search_dir,alpha))
            alpha *= 1/2
            if alpha < 10^-6
                accept = false
                break
            end
        end
        return alpha,accept
    end

    #---------------------#---------------------
    #BFGS Hessian approximation.
    function BFGS(current_p_vars,grad)
        l = length(prev_p_vars)
        if l < 2
            return [1 0; 0 1]
        end    
        Hk = prev_H   
        dk = grad-prev_grad
        pk = current_p_vars-prev_p_vars[l-1]

        return Hk - (Hk*dk)*(Hk*dk)'/(dk'*Hk*dk)+(pk*pk')/(pk'*dk)
    end

    #---------------------#---------------------
    #Calculate step (direction and length).
    function calc_step(current_p_vars,grad)
        # hessian_approx = (1/2).*grad*grad'      <--- old hessian approx
        hessian_approx = BFGS(current_p_vars,grad)
        prev_H = hessian_approx
        search_dir = -inv(hessian_approx)*grad #Newton search direction

        alpha,accepted = armijo(current_p_vars,grad,search_dir)

        if !accepted
            search_dir = -grad
            alpha = steepest_descent_step(current_p_vars,search_dir)
        end

        return alpha*search_dir
    end

    #---------------------#---------------------
    #Termination criteria.
    function term(grad,step)
        i = length(prev_costs)
        if i<2
            return false
        end

        C1 = norm(grad) <= epsilon[1]*(1+abs(prev_costs[i]))
        C2 = prev_costs[i-1] - prev_costs[i] <= epsilon[2]*(1+abs(prev_costs[i]))
        C3 = norm(step) <= epsilon[3]*(1+norm(prev_p_vars[i]))

        return (C1 && C2) || (C2 && C3) || (C1 && C3)
    end

    #---------------------#---------------------
    #Estimation loop.
    while true
        #Add current parameters/cost to "memory".
        push!(prev_p_vars,p_var)
        push!(prev_costs,cost_fun_of_p(p_var))

        steady_state_sol = Steady_state_solver(p_const,p_var,pre_shift)

        sensitivity_sol = sensitivity_solver(steady_state_sol,post_shift,tvals,p_const,p_var)

        #Calculate gradient and step, then take the step.
        gradient = calc_grad_sens(sensitivity_sol)
        step = calc_step(p_var,gradient)
        p_var += step
        prev_grad = gradient

        #Termination.
        if term(gradient,step)
            break
        end
    end
    
    #Add the found parameters and cost to "memory", then return "memory" vectors.
    push!(prev_p_vars,p_var)
    push!(prev_costs,cost_fun_of_p(p_var))
    return prev_p_vars, prev_costs
end

"""
    initiate_constants()

    Initiate the constants epsilon, mu and nu for termination criteria and step length calculations.
"""
function initiate_constants()
    epsilon = [10^-6,10^-5,10^-3] #termination criteria
    mu = 0.01 #Wolfe 1 (aka Armijo)
    nu = 0.9 #Wolfe 2
    return epsilon,mu,nu
end

"""
    cost_fun(prediction)

    Calculate cost of prediction.
"""
function cost_fun(pred)
    s = 0
    for t in 1:length(pred)
        for y in 1:length(pred[1])
            s += (pred[t][y]-data[t][y])^2
        end
    end
    return s
end

"""
    calc_grad_sens(sensitivity_sol)

    Calculate gradient using sensitivity analysis.
"""
function calc_grad_sens(sens_sol)
    sens_mat = extract_local_sensitivities(sol)

    #y(p,t) = #sens_sol.u[t][y]
    #dy_dp = sens_mat[2][p][y,t]
    n_ps = size(sens_mat[2])[0]
    n_ts,n_ys = size(sens_sol.u)

    grad = zeros(n_ps)
    if sol.retcode != :Success
        print("Instability... Setting gradient to 0.\n") #ooga booga
        return grad
    end

    for p = 1:n_ps
        for y = 1:n_ys
            for t = 1:n_ts
                grad[p] += (sens_sol.u[t][y]-data[t][y])*sens_mat[2][p][y,t]
            end
        end
    end
    return grad.*2
end

"""
    calc_grad_AD(...(TBD)...)

    Calculate gradient using automatic differentiation.
"""
function calc_grad_AD(p_var, p_const)
    function dMK_dp_Glucose_addition_Mig1(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 20) #ODE for MIG1 on row 20 in ODE system

            MK_for_diff = sum((data_for_ode - (log.(10,u_approx./(1e-5 .+ 1.0 .- u_approx)))).^2)
            return MK_for_diff
        end
    
        data_for_ode = data__Mig1_glucose_relief
        timelist_for_ode = t_Mig1_glucose_relief
        tspan = (0.0, 20.0) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_var, p_const, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) 
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Glucose_addition_Sch9(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 6)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_sch9Delta_cAMP, 0.052, 1.227)
        timelist_for_ode = t_sch9Delta_cAMP
        tspan = (0.0, 3.5) # [min]
        p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) 
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]
        p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 1.0

        return dMK_dp
    end

    function dMK_dp_Glucose_addition_cAMP(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 6) 

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_cAMP, 0.052, 1.227)
        timelist_for_ode = t_cAMP
        tspan = (0.0, 3.5) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Glucose_addition_Sch9_p(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12) 

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Sch9_glucose_relief)
        timelist_for_ode = t_Sch9_glucose_relief
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Glucose_starvation_Snf1_p(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 10)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Snf1)
        timelist_for_ode = t_Snf1
        tspan = (0.0, 62) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Glucose_starvation_Sch9_p(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Sch9_glucose_starve)
        timelist_for_ode = t_Sch9_glucose_starve
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Sch9P_glutamine_L(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Sch9P_glutamine_L, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_L
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Sch9P_glutamine_H(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Sch9P_glutamine_H, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_H
        tspan = (0.0, 32) # [min]
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]

        return dMK_dp
    end

    function dMK_dp_Glutamine_addition_Sch9_gtr1Delta(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 12)

            MK_for_diff = sum((data_for_ode - u_approx).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Sch9_gtr1Delta, 4.82, 52.27)
        timelist_for_ode = t_Sch9_gtr1Delta
        tspan = (0.0, 32) # [min]
        w_torc_ego_true = p_var[Get_index(p_var_lookup_table, "w_torc_ego")]
        w_torc_egoin_true = p_var[Get_index(p_var_lookup_table, "w_torc_egoin")]

        # Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
        p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
        p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
        p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]
        p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 1.0
        p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego_true
        p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin_true

        return dMK_dp
    end

    function dMK_dp_Rapamycin_treatment(p_var,p_const)
        function calc_FWD_grad(p_values)
            p_values = exp.(p_values)
            sol = ODE_solver_FWDgrad(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_values, timelist_for_ode, prob)
            u_approx = getindex.(sol.u, 25)

            MK_for_diff = sum((data_for_ode - (u_approx./(1e-3 .+ u0_SS[25]))).^2)
            return MK_for_diff
        end
    
        data_for_ode = Minmaxnorm(data_Rib_rap)
        timelist_for_ode = t_Rib_rap
        tspan = (0.0, 92.0) # [min]
        p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0
        p_values = vcat(p_var,p_const)
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
        prob = ODEProblem(ODE_sys, u0_SS, tspan, p_values)
        dMK_dp = ForwardDiff.gradient(calc_FWD_grad, log.(prob.p))[1:length(last.(p_var))]
        p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 1.0

        return dMK_dp
    end

    dMK_dp_Glucose_addition_Mig1_solution = dMK_dp_Glucose_addition_Mig1(p_var,p_const)
    dMK_dp_Glucose_addition_Sch9_solution = dMK_dp_Glucose_addition_Sch9(p_var,p_const)
    dMK_dp_Glucose_addition_cAMP_solution = dMK_dp_Glucose_addition_cAMP(p_var,p_const)
    dMK_dp_Glucose_addition_Sch9_p_solution = dMK_dp_Glucose_addition_Sch9_p(p_var,p_const)
    dMK_dp_Glucose_starvation_Snf1_p_solution = dMK_dp_Glucose_starvation_Snf1_p(p_var,p_const)
    dMK_dp_Glucose_starvation_Sch9_p_solution = dMK_dp_Glucose_starvation_Sch9_p(p_var,p_const)
    dMK_dp_Sch9P_glutamine_L_solution = dMK_dp_Sch9P_glutamine_L(p_var,p_const)
    dMK_dp_Sch9P_glutamine_H_solution = dMK_dp_Sch9P_glutamine_H(p_var,p_const)
    dMK_dp_Glutamine_addition_Sch9_gtr1Delta_solution = dMK_dp_Glutamine_addition_Sch9_gtr1Delta(p_var,p_const)
    dMK_dp_Rapamycin_treatment_solution = dMK_dp_Rapamycin_treatment(p_var,p_const)

    True_grad = dMK_dp_Glucose_addition_Mig1_solution + dMK_dp_Glucose_addition_Sch9_solution + dMK_dp_Glucose_addition_cAMP_solution + dMK_dp_Glucose_addition_Sch9_p_solution + dMK_dp_Glucose_starvation_Snf1_p_solution + dMK_dp_Glucose_starvation_Sch9_p_solution + dMK_dp_Sch9P_glutamine_L_solution + dMK_dp_Sch9P_glutamine_H_solution + dMK_dp_Glutamine_addition_Sch9_gtr1Delta_solution + dMK_dp_Rapamycin_treatment_solution

    return True_grad
end