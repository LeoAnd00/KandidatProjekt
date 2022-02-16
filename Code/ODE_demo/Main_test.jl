# Löser ODEer i ODE_functions, samma index som där med [1]=x, [2]=y, [3]=z
using DifferentialEquations
include("ODE_functions_test.jl"), ("Get_parmeter_values_test.jl")
u0=[1.0, 0.0, 0.0]
tspan=(0.0, 100.0)
prob=ODEProblem(ODE_test, u0, tspan)
sol=solve(prob)