using DifferentialEquations
using LinearAlgebra

"""
    est_param(ODE_system, tvals, initial_values, start_guess[, n_parameters])

    Estimate parameters for ODE_system.

    If n_parameters is unspecified, estimation is done for all parameters in start_guess.
"""
function est_param(ODE_sys, tvals, init_vals, start_guess, n_ps = length(start_guess))

    param = start_guess[1:n_ps]
    tspan = [first(tvals),last(tvals)]

    #Information updated in each iteration.
    prev_params = [] #list of all previous parameter values, "memory"
    prev_costs = [] #list of corresponding costs, "memory"
    prev_H = [1 0; 0 1] #last Hessian approximation
    prev_grad = [0,0] #last gradient

    #Constants: epsilon used in termination criteria and mu, nu used in step length calculations.
    epsilon,mu,nu = initiate_constants()

    #---------------------#---------------------
    #Cost as function of parameters. Uses cost_fun(pred) below.
    function cost_fun_of_param(params)
        prob = ODEProblem(ODE_sys,init_vals,tspan,params)
        sol = solve(prob,Rodas5(),saveat=tvals,verbose = false)
        if sol.retcode != :Success
            return Inf
        end
        return cost_fun(sol.u)
    end

    #---------------------#---------------------
    #Conditions used when calculating the step length alpha.

    #General condition on alpha: positive parameters.
    function GC(current_params,search_dir,alpha)
        for i in 1:n_ps
            if current_params[i]+alpha*search_dir[i] < 0
                return false
            end
        end
        return true
    end

    #Condition on alpha used for steepest descent: new cost lower than previous cost.
    function SDC(current_params,search_dir,alpha)
        return cost_fun_of_param(current_params+alpha*search_dir) <= cost_fun_of_param(current_params)
    end

    #Wolfe condition 1 (aka Armijo condition).
    function WC1(current_params,grad,search_dir,alpha)
        return cost_fun_of_param(current_params+alpha*search_dir)-cost_fun_of_param(current_params) <= mu*alpha*search_dir'*grad
    end

    #Wolfe condition 2 (strong form).
    function WC2(current_params,grad,search_dir,alpha)
        return abs(search_dir'*calc_grad(exp_ODE_sys,init_vals,tspan,current_params+alpha*search_dir,tvals)) >= abs(nu*search_dir'*grad)
    end

    #---------------------#---------------------
    #Functions to calculate step length (alpha).

    #For steepest descent.
    function steepest_descent_step(current_params,search_dir)
        alpha = 1
        while !(SDC(current_params,search_dir,alpha) && GC(current_params,search_dir,alpha))
            alpha *= 1/2
        end
        return alpha
    end

    #Armijo step length rule (using Wolfe conditions).
    function armijo(current_params,grad,search_dir)
        alpha = 1
        accept = true
        while !(WC1(current_params,grad,search_dir,alpha) && WC2(current_params,grad,search_dir,alpha) && GC(current_params,search_dir,alpha))
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
    function BFGS(current_params,grad)
        l = length(prev_params)
        if l < 2
            return [1 0; 0 1]
        end    
        Hk = prev_H   
        dk = grad-prev_grad
        pk = current_params-prev_params[l-1]

        return Hk - (Hk*dk)*(Hk*dk)'/(dk'*Hk*dk)+(pk*pk')/(pk'*dk)
    end

    #---------------------#---------------------
    #Calculate step (direction and length).
    function calc_step(current_params,grad)
        # hessian_approx = (1/2).*grad*grad'      <--- old hessian approx
        hessian_approx = BFGS(current_params,grad)
        prev_H = hessian_approx
        search_dir = -inv(hessian_approx)*grad #Newton search direction

        alpha,accepted = armijo(current_params,grad,search_dir)

        if !accepted
            search_dir = -grad
            alpha = steepest_descent_step(current_params,search_dir)
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
        C3 = norm(step) <= epsilon[3]*(1+norm(prev_params[i]))

        return (C1 && C2) || (C2 && C3) || (C1 && C3)
    end

    #---------------------#---------------------
    #Estimation loop.
    while true
        #Add current parameters/cost to "memory".
        push!(prev_params,param)
        push!(prev_costs,cost_fun_of_param(param))

        #Calculate gradient and step, then take the step.
        gradient = calc_grad_sens(exp_ODE_sys,init_vals,tspan,param,tvals)
        step = calc_step(param,gradient)
        param += step
        prev_grad = gradient

        #Termination.
        if term(gradient,step)
            break
        end
    end
    
    #Add the found parameters and cost to "memory", then return "memory" vectors.
    push!(prev_params,param)
    push!(prev_costs,cost_fun_of_param(param))
    return prev_params, prev_costs
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
    calc_grad_sens(...(TBD)...)

    Calculate gradient using sensitivity analysis. (Trasig pga n_ts, etc.)
"""
function calc_grad_sens(ODE_exp,init_vals,tspan,param,tvals)
    init_vals = [init_vals; ones(4)]                                #ooga booga
    prob = ODEProblem(ODE_exp,init_vals,tspan,param)
    sol = solve(prob,Rodas5(),saveat=tvals,verbose = false)

    grad = zeros(n_ps)

    if sol.retcode != :Success
        print("Instability... Setting gradient to 0.\n")            #ooga booga
        return grad
    end

    for p = 1:n_ps
        for y = 1:n_ys
            for t = 1:n_ts
                grad[p] += (sol.u[t][y]-data[t][y])*sol.u[t][y*2+p] #ooga booga
            end
        end
    end

    return grad.*2
end