# Definierar och löser ODE-problemet

using DifferentialEquations
using ModelingToolkit
using Plots

include("ODE_functions.jl") 
u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 
tspan = (0.0, 70.0) # [min]

include("parameter_values.jl")
prob = ODEProblem(ODE_sys, u0, tspan, p) # Definierar ODE-problemet

sol = solve(prob) # Löser ODE-problemet
plot(sol, legend = false) # Plottar lösningen, utan legend