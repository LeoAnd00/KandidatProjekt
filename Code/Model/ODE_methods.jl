# Innehåller de funktioner och verktyg som är nödvändiga för att lösa steady-state och ODE-systemet 
# Innehåller även ett verktyg för att normalisera data och ett verktyg för att hitta index i lookup_table

using DifferentialEquations
using DiffEqSensitivity

include("ODE_functions.jl") 

# Steady-state_solver löser ut steady-state för ODE-systemet vid givna parametervärden och näringsvärden
# INPUTS: model_inputs ges som en vektor av par, enligt (name => value, ...), där alla inputs som inte är noll behöver vara med. 
# OUTPUT: En vektor av steady-state värden 

function Steady_state_solver(p_const, p_var, model_inputs)
    
    p = vcat(p_const, p_var)

    u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

    for i in range(1, length(model_inputs))
        u0[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end

    SS_prob = SteadyStateProblem(ODE_sys, u0, p) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem

    SS_sol = solve(SS_prob, DynamicSS(Rodas5(), abstol=1e-9, reltol=1e-9)) # Löser steady-state problemet
    return SS_sol.u
end

function Steady_state_solver(p_conc, model_inputs)
    
    u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

    for i in range(1, length(model_inputs))
        u0[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end
    SS_prob = SteadyStateProblem(ODE_sys, u0, p_conc) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem

    SS_sol = solve(SS_prob, DynamicSS(Rodas5(), abstol=1e-9, reltol=1e-9)) # Löser steady-state problemet
    return SS_sol.u
end

function ODE_solver(u0_SS, model_inputs, tspan, p_const, p_var)

    p = vcat(p_const, p_var) 

    for i in range(1, length(model_inputs))
        u0_SS[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end

    prob = ODEProblem(ODE_sys, u0_SS, tspan, p)
    return solve(prob, Rodas5(), abstol=1e-9, reltol=1e-9)
end

function ODE_solver(u0_SS, model_inputs, tspan, p_conc)

    for i in range(1, length(model_inputs))
        u0_SS[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end

    prob = ODEProblem(ODE_sys, u0_SS, tspan, p_conc)
    return solve(prob, Rodas5(), abstol=1e-9, reltol=1e-9)
end

function sensitivity_solver(u0_SS, model_inputs, tvals, p_const, p_var)

    p = vcat(p_const, p_var)

    for i in range(1, length(model_inputs))
        u0_SS[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end

    prob = ODEForwardSensitivityProblem(ODE_sys, u0_SS, [first(tvals),last(tvals)], p, sensealg=ForwardSensitivity(autodiff=false))

    return solve(prob,Rodas5(autodiff = false),saveat=tvals, sensealg=ForwardSensitivity(autodiff=false))
end

function sensitivity_solver(u0_SS, model_inputs, tspan, p_conc)

    for i in range(1, length(model_inputs))
        u0_SS[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end

    prob = ODEForwardSensitivityProblem(ODE_sys, u0_SS, tspan, p_conc)

    return solve(prob, Rodas5())
end

# Minmaxnorm normaliserar värdena i input list.
# Finns i två varianter: En där min och max hämtas från listan, och en där de matas in manuellt

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

function Steady_state_solver_FWD(p_const_FWD, p_var_FWD, model_inputs)
    
    p = vcat(p_var_FWD, p_const_FWD)

    u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

    for i in range(1, length(model_inputs))
        u0[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    SS_prob = SteadyStateProblem(ODE_sys, u0, p) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem

    SS_sol = solve(SS_prob, DynamicSS(Rodas4P(), abstol=1e-3, reltol=1e-3), abstol=1e-8, reltol=1e-8)#,  maxiters = 1e7) # Löser steady-state problemet

    #println(" result_steady_state: ", SS_sol.retcode)

    steady_state_solution = SS_sol.u

    # Checks derivatives
    if SS_sol.retcode != :Success
        du = zeros(length(SS_prob.u0))
        SS_prob.f(du, steady_state_solution, SS_prob.p, 1.0)
        #println(" all: ",du)
        for i in 1:length(du)
            if i ∉ [2,3,6,8,11]
                if abs(du[i]) > 0.03
                    println("i that gave to big gradient after steady state")
                    println("i = $i")
                    steady_state_solution = "Not working"# Will be flagged as an error in the next step of optimization
                    break
                end
            end
        end
    end

    return steady_state_solution
end

function ODE_solver_FWDgrad(u0_SS, model_inputs, tspan, p_values, timelist_for_ode, prob)

    prev_model_inputs = (Carbon => u0_SS[end-2], ATP => u0_SS[end-1], Glutamine_ext => u0_SS[end])

    for i in range(1, length(model_inputs))
        u0_SS[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    nprob = remake(prob; p = p_values)
    nprob = remake(nprob; u0_SS = u0_SS)
    nprob = remake(nprob; tspan = tspan)
    sol = solve(nprob,Rodas4P(),abstol=1e-8,reltol=1e-8, saveat = timelist_for_ode)

    for i in range(1, length(prev_model_inputs))
        u0_SS[findfirst(x->x==string(first.(prev_model_inputs)[i]), u_lookup_table)] = last.(prev_model_inputs)[i]
    end

    return sol#, maxiters = 1e5,dtmin = 1e-5)
end

function ODE_solver_FWD(u0_SS, model_inputs, tspan, p_const_FWD, p_var_FWD, timelist_for_ode)

    p = vcat(p_var_FWD, p_const_FWD) 

    prev_model_inputs = (Carbon => u0_SS[end-2], ATP => u0_SS[end-1], Glutamine_ext => u0_SS[end])

    for i in range(1, length(model_inputs))
        u0_SS[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    prob = ODEProblem(ODE_sys, u0_SS, tspan, p)

    for i in range(1, length(prev_model_inputs))
        u0_SS[findfirst(x->x==string(first.(prev_model_inputs)[i]), u_lookup_table)] = last.(prev_model_inputs)[i]
    end

    return solve(prob,Rodas4P(),abstol=1e-8,reltol=1e-8, saveat = timelist_for_ode)#, maxiters = 1e5,dtmin = 1e-5)
end

nutrient_shifts = [
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0) => (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0) => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0) => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), 
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0) => (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)
]

nutrient_shifts_lookup_table = ["Glucose starvation", "Glucose addition", "High glutamine", "Nitrogen starvation"] # Shiftningarnas position

index_glucose_starvation = Get_index(nutrient_shifts_lookup_table, "Glucose starvation")
index_glucose_addition = Get_index(nutrient_shifts_lookup_table, "Glucose addition")
index_high_glutamine = Get_index(nutrient_shifts_lookup_table, "High glutamine")
index_nitrogen_starvation = Get_index(nutrient_shifts_lookup_table, "Nitrogen starvation")