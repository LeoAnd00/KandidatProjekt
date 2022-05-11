# Starts by making a latin hypercube sample with 10000 samples on 
# logarithmic scale between -3 and 3, with x working ones. Then these are 
# optimized.
include("../Latin_Hypercube_Sampling/latin_hypercube_sampling.jl")
main()
include("optimize_all_p_0_values.jl")