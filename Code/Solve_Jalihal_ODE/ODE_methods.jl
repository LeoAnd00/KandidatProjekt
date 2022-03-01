using DifferentialEquations
using ModelingToolkit

include("ODE_functions.jl") 

function ODE_solver(u0_SS, tspan, p)
    prob = ODEProblem(ODE_sys, u0_SS, tspan, p)
    return solve(prob)
end

function Steady_state_solver(p_preshift)
    
    u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

    SS_prob = SteadyStateProblem(ODE_sys, u0, p_preshift) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem

    SS_sol = solve(SS_prob) # Löser steady-state problemet
    return SS_sol.u
end