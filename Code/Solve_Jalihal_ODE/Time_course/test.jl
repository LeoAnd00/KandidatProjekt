model_inputs = (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)

u0 = zeros(length(eqs)) # Initialkoncentrationer för variabler är noll 

for i in range(1, length(model_inputs))
    u0[findfirst(x->x==string(first.(model_inputs)[i]), u_lookup_table)] = last.(model_inputs)[i]
end
SS_prob = SteadyStateProblem(ODE_sys, u0, p_preshift) # Definierar steady-state problemet, vilket är ett simpelt ODE-problem
