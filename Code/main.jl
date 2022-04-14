using LatinHypercubeSampling
include("Sensitivity_Analysis/sensitivity_analysis.jl")

n_ps = 81
scale = repeat([(-3,3)],n_ps)
p0 = scaleLHC(randomLHC(100,n_ps),scale)

sensitivities = get_sensitivity_matrix(exp.(p0[1,:]))
