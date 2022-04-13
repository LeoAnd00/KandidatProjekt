using LatinHypercubeSampling
include("Parameter_Estimation/param_estim.jl")
include("Solve_Jalihal_ODE/Time_course/ODE_methods.jl")

#p0 via latin hypercube på logskala (log 10) från -3 till 3
n_ps = 2 #ska vara antal p_var
p0 = scaleLHC(randomLHC(100,n_ps),[(-3,3),(-3,3)])