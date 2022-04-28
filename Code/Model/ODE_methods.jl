# Innehåller de funktioner och verktyg som är nödvändiga för att lösa steady-state och ODE-systemet 
# Innehåller även ett verktyg för att normalisera data och ett verktyg för att hitta index i lookup_table

using DifferentialEquations
using DiffEqSensitivity
using Documenter

include("ODE_functions.jl") 

"""
    Steady_state_solver(p_const, p_var, model_inputs)
    Steady_state_solver(p_conc, model_inputs)

Compute the steady state for the parameters and nutrient conditions in `model_input`

# Examples
```julia-repl
julia> Steady_state_solver(p_conc, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))
```
"""
function Steady_state_solver(p_const, p_var, model_inputs)
    
    p_cat = vcat(p_const, p_var) # Concatenates the component vectors in the form required as input by ModelingToolkit

    u0 = zeros(length(eqs))     # Initial concentrations of are zero

    # Implements the nutrient condition values as initial guesses, values will be kept during calculations
    for i in range(1, length(model_inputs))
        u0[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end

    SS_prob = SteadyStateProblem(ODE_sys, u0, p_cat)

    SS_sol = solve(SS_prob, DynamicSS(Rodas5(), abstol=1e-9, reltol=1e-9))
    return SS_sol.u
end

function Steady_state_solver(p_conc, model_inputs)
    
    u0 = zeros(length(eqs)) 

    for i in range(1, length(model_inputs))
        u0[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
    end
    SS_prob = SteadyStateProblem(ODE_sys, u0, p_conc) 

    SS_sol = solve(SS_prob, DynamicSS(Rodas5(), abstol=1e-9, reltol=1e-9)) 
    return SS_sol.u
end

"""
    ODE_solver(u0_SS, model_inputs, tspan, p_const, p_var)
    ODE_solver(u0_SS, model_inputs, tspan, p_conc)

Compute the solution of the ODE system with nutrient conditions of `model_inputs`, and the initial conditions in `u0_SS`.
Unit of time in `tspan` is minutes. 

# Examples
```julia-repl
julia> ODE_solver(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_conc)
```
"""
function ODE_solver(u0_SS, model_inputs, tspan, p_const, p_var)

    p = vcat(p_const, p_var) 

    # Implements the nutrient conditions as initial conditions, values will be kept during the calculations
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

"""

"""
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

"""
    Minmaxnorm(input_list)
    Minmaxnorm(input_list, min, max)

Normalize values in `input_list` according to a minimum and maximum value.
For an element in the list, normalization occurs according to (element-min)/(max-min).
The min/max values is either retrieved automatically from `input_list` or as function inputs.

# Examples
```julia-repl
julia> Minmaxnorm([4, 2, 4, 1])
4-element Vector{Any}:
 1.0
 0.3333333333333333
 1.0
 0.0
```
```julia-repl
julia> Minmaxnorm([4, 2, 4, 1], 0.5, 6)
4-element Vector{Any}:
0.6363636363636364
0.2727272727272727
0.6363636363636364
0.09090909090909091
```
"""
function Minmaxnorm(input_list)
    min = minimum(input_list)
    max = maximum(input_list)

    output_list = []
    for element in input_list
        append!(output_list, (element - min)/(max-min))     # Fills the output_list with normalized values
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

"""
    Get_index(list_lookup_table, key)

Retrieve the index of an element in a list using the corresponding `list_lookup_table` and the `key` in the form of a string.
"""
function Get_index(list_lookup_table, key)
    return findfirst(x->x==key, list_lookup_table)    # Finds the first instance of the key and returns the index 
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

"""

"""
function Steady_state_solver_FWD(p_const_FWD, p_var_FWD, model_inputs)
    
    p = vcat(p_var_FWD, p_const_FWD)

    u0 = zeros(length(eqs)) 

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

"""

"""
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

"""

"""
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