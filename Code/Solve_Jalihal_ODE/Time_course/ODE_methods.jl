using DifferentialEquations
using ModelingToolkit

include("ODE_functions.jl") 

function Steady_state_solver(p_preshift, model_inputs)
    
    u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

    for i in range(1, length(model_inputs))
        u0[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end
    SS_prob = SteadyStateProblem(ODE_sys, u0, p_preshift) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem

    SS_sol = solve(SS_prob, DynamicSS(Rodas5(), abstol=1e-12, reltol=1e-12)) # Löser steady-state problemet
    return SS_sol.u
end

function ODE_solver(u0_SS, model_inputs, tspan, p)

    for i in range(1, length(model_inputs))
        u0_SS[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
    end

    prob = ODEProblem(ODE_sys, u0_SS, tspan, p)
    # return solve(prob, Rodas5(), abstol=1e-12, reltol=1e-12)
    return solve(prob)
end

