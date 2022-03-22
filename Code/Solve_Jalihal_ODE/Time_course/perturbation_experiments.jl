#
include("perturbation_methods.jl")

# Steady_state_pertubation_solver()

# Snf1
# Pre-shift => Glutamine_ext = 1,Carbon=1, ATP=1
#p[findfirst(x->x=="Sak_T", p_lookup_table)] = Sak_T => 0.0

include("ODE_methods.jl")
# u0_SS_Snf1 = Steady_state_solver(p, (Carbon => 1.0, Glutamine_ext => 1.0, ATP => 1.0)) # Returnerar steady state för parametrarna p

# # POST => Carbon = 0.0, ATP = 0.0

# u0_SS_Snf1_2 = Steady_state_solver(p, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0))

# println("Pre-shift: ", u0_SS_Snf1[10], ", post-shift: ", u0_SS_Snf1_2[10])

# Steady_state_pertubation_solver((Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), Snf1)

Steady_state_pertubation_solver_mut([Sch9_T => 0.0], (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), Gis1)