using LatinHypercubeSampling
using LinearAlgebra
using DifferentialEquations
using DiffEqSensitivity
using ModelingToolkit
using CSV
using DataFrames
include("../Model/ODE_functions.jl")
include("../Model/parameter_values.jl")
include("../Parameter_Estimation/Optimization.jl")
######################################### LatinHypercubeSampling ################################################################
"""
    Check_if_p_works(p_var, p_const, scaled_plan)

    Checks if the parameter sets from the latin hypercube works and saves them into a CSV file.
"""
function Check_if_p_works(p_var, p_const, scaled_plan)

    N = 0
    df = DataFrame()

    for item in 1:length(scaled_plan[:,1])

        print(" Iteration: ",item)

        p_values_process = exp10.(scaled_plan[item,:])
        #p_values_process = scaled_plan[item,:]
        p_values_dual_temp = [Pair(first.(p_var)[1], p_values_process'[1])]
        for i in range(2, length(last.(p_var)))
            B = Pair(first.(p_var)[i], p_values_process'[i])
            p_values_dual_temp = vcat(p_values_dual_temp,B)
        end
        p_var_new = p_values_dual_temp

        # Check if the cost can be calculated
        MK = 0
        try
            MK = calc_cost(p_var_new,p_const)
        catch
            continue
        end

        N += 1
        # Saves parameters that works in a CSV file
        if N == 1
            df = DataFrame(Parameter = vcat(first.(p_var_new),Num("Cost")),
                        p_0_value_1 = vcat(last.(p_var_new),MK)) 
            CSV.write("Intermediate/p_0_values.csv",df)
        else
            df[!,Symbol("p_0_value_$N")] = vcat(last.(p_var_new),MK)
            CSV.write("Intermediate/p_0_values.csv",df)
        end

    end

    return N
end

"""
    main()

    Here you define different settings for the LHS, like number of samples, and then it generates a LHS and runs the Check_if_p_works(p_var, p_const, scaled_plan)
"""
function main()
    include("../Model/parameter_values.jl")

    plan = randomLHC(10000,81) #samples, dimensions

    #plan_scale = [(0.95*last.(p_var)[1],1.05*last.(p_var)[1])]
    plan_scale = [(-3.0,3.0)]
    for i in 1:80
        #push!(plan_scale,(0.95*last.(p_var)[i+1],1.05*last.(p_var)[i+1]))
        push!(plan_scale,(-3.0,3.0))
    end

    scaled_plan = scaleLHC(plan,plan_scale)

    N = Check_if_p_works(p_var, p_const, scaled_plan)
    println(" Number of working parameter sets: ", N) 

end

main()

################################################################################################################################
