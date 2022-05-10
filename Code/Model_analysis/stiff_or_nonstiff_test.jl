using DifferentialEquations
include("../Model/ODE_functions.jl")
include("../Model/ODE_methods.jl")
include("../Model/parameter_values.jl")

u0_SS = Steady_state_solver(p_conc, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0))

model_inputs = (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)

for i in range(1, length(model_inputs))
    u0_SS[Get_index(u_lookup_table, string(first.(model_inputs)[i]))] = last.(model_inputs)[i]
end

prob = ODEProblem(ODE_sys, u0, (0.0, 30.0), p_conc)

# Compares solve times of stiff and non-stiff solver
println("Stiff solver:")
@time sol_stiff = solve(prob, Kvaerno4(), abstol=1e-8, reltol=1e-8)
println("Non-stiff solver:")
@time sol_nonstiff = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)
println("RodasS (stiff-aware):")
@time sol_Rodas5 = solve(prob, Rodas5(), abstol=1e-8, reltol=1e-8)
println("Rodas4P (stiff-aware):")
@time sol_Rodas4P = solve(prob, Rodas4P(), abstol=1e-8, reltol=1e-8)

# The non-stiff solver ends up being faster than the stiff solver => The ODE problem is stiff
# In the end, Rodas4P() ends up being the fastest solver