using LatinHypercubeSampling
using LinearAlgebra
using DifferentialEquations
using DiffEqSensitivity
using ModelingToolkit
using CSV
using DataFrames
include("../Model/ODE_functions.jl")
include("../Model/parameter_values.jl")
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

        #p_values_process = exp10.(scaled_plan[item,:])
        p_values_process = scaled_plan[item,:]
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

    #df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
    #print(df)

    return N
end

"""
    Check_if_p_MK_under_50(x)

    Makes CSV with working parameter sets with MK under x
"""
function Check_if_p_MK_under_x(x)

    N = 0
    df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
    for i in 1:(length(df[1,:])-1)
        df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
        MK_df = df[:, Symbol("p_0_value_$i")][end]
        MK_that_works = df[:, Symbol("p_0_value_$i")][1:end]
        if MK_df > x
            continue
        else
            N += 1
            if N == 1
                df = DataFrame(Parameter = vcat(first.(p_var),Num("Cost")),
                            p_0_value_1 = MK_that_works)
                CSV.write("Intermediate/p_0_values_MK_under_x.csv",df)
            else
                df = CSV.read("Intermediate/p_0_values_MK_under_x.csv", DataFrame)
                df[!,Symbol("p_0_value_$N")] = MK_that_works
                CSV.write("Intermediate/p_0_values_MK_under_x.csv",df)
            end
        end
    end
    return N
end

function main()
    include("../Model/parameter_values.jl")

    plan = randomLHC(100,81) #samples, dimensions

    #plan_scale = [(0.7*last.(p_var)[1],1.5*last.(p_var)[1])]
    plan_scale = [(0.95*last.(p_var)[1],1.05*last.(p_var)[1])]
    #plan_scale = [(-3.0,3.0)]
    #plan_scale = [(-3.0,2.2)]
    #plan_scale = [(0.001,120)]
    for i in 1:80
        push!(plan_scale,(0.95*last.(p_var)[i+1],1.05*last.(p_var)[i+1]))
        #push!(plan_scale,(0.7*last.(p_var)[i+1],1.5*last.(p_var)[i+1]))
        #push!(plan_scale,(-3.0,3.0))
        #push!(plan_scale,(-3.0,2.2))
        #push!(plan_scale,(0.001,120))
    end

    scaled_plan = scaleLHC(plan,plan_scale)

    N = Check_if_p_works(p_var, p_const, scaled_plan)
    println(" Number of working parameter sets: ", N) 

    constraint_under_x = false
    if constraint_under_x == true
        x = 50
        N2 = Check_if_p_MK_under_x(x)
        println(" Number of working parameter sets under $x MK: ", N2) 
    end
end

main()

################################################################################################################################
