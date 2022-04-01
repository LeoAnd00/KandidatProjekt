function Steady_state_solver(p_const, p_var, model_inputs)
    
    p = vcat(p_var, p_const)

    u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

    for i in range(1, length(model_inputs))
        u0[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    SS_prob = SteadyStateProblem(ODE_sys, u0, p) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem

    SS_sol = solve(SS_prob, DynamicSS(Rodas5(), abstol=1e-12, reltol=1e-12)) # Löser steady-state problemet
    return SS_sol.u
end

function ODE_solver(u0_SS, model_inputs, tspan, p_const, p_var, timelist_for_ode)

    p = vcat(p_var, p_const) 

    for i in range(1, length(model_inputs))
        u0_SS[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    prob = ODEProblem(ODE_sys, u0_SS, tspan, p)

    return solve(prob,Rodas5(),abstol=1e-8,reltol=1e-8, saveat = timelist_for_ode)#, maxiters = 1e5,dtmin = 1e-5)
end

function Minmaxnorm(input_list)
    min = minimum(input_list)
    max = maximum(input_list)

    output_list = []

    for element in input_list
        append!(output_list, (element - min)/(max-min))
    end

    return output_list
end

function Minmaxnorm(input_list, min, max)

    output_list = []

    for element in input_list
        append!(output_list, (element - min)/(max-min))
    end

    return output_list
end

function Get_index(input_list, key)
    return findfirst(x->x==key, input_list)     
end

function ODE_solver_grad(u0_SS, model_inputs, tspan, p_values, timelist_for_ode, prob)
            
    for i in range(1, length(model_inputs))
        u0_SS[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    nprob = remake(prob; p = p_values)
    nprob = remake(nprob; u0_SS = u0_SS)
    nprob = remake(nprob; tspan = tspan)

    return solve(nprob,Rodas5(),abstol=1e-8,reltol=1e-8, saveat = timelist_for_ode)#, maxiters = 1e5,dtmin = 1e-5)
end