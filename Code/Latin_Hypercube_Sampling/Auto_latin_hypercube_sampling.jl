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
    main()

    Finds the best parametervector in a list called p_0_values.csv and then generates latin hypercube samples around it
    (user defines interval), takes the best result and repeats until the user stops the iteration. If none of the costs were 
    below the previously best one it will try to increase the number of samples by a choosen amount and generate a new latin hypercube.
"""
function main()
    include("../Model/parameter_values.jl")
    MK_current_best = 1e100
    MK_list = []
    p_0_var = p_var
    N = 0
    Number_of_samples = 70
    while true
        N += 1
        p_var = p_0_var
        if N > 1
            plan = randomLHC(Number_of_samples,81) #samples, dimensions

            plan_scale = [(0.97*last.(p_var)[1],1.03*last.(p_var)[1])]
            #plan_scale = [(-3.0,3.0)]
            for i in 1:80
                push!(plan_scale,(0.97*last.(p_var)[i+1],1.03*last.(p_var)[i+1]))
                #push!(plan_scale,(-3.0,3.0))
            end

            scaled_plan = scaleLHC(plan,plan_scale)

            N_check = Check_if_p_works(p_var, p_const, scaled_plan)
            println(" Number of working parameter sets: ", N_check)  
        end

        MK_df = 1e100
        Best_i = 0
        df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
        for i in 1:(length(df[1,:])-1)
            df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
            MK_df_new = df[:, Symbol("p_0_value_$i")][end]
            if MK_df_new < MK_df && MK_df_new < MK_current_best
                MK_df = MK_df_new
                Best_i = i
            end
        end
        println("Best_i:",Best_i)
        if Best_i == 0
            Number_of_samples = Number_of_samples + 20
            continue
        end
        try
            df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
            p_names_var = first.(p_var)
            p_0_var = df[:, Symbol("p_0_value_$Best_i")][1:end-1]

            p_var_temp = [Pair(p_names_var[1], p_0_var'[1])]
            for item in range(2, length(p_0_var))
                B = Pair(p_names_var[item], p_0_var'[item])
                p_var_temp = vcat(p_var_temp,B)
            end
            p_0_var = p_var_temp
            MK_current_best = MK_df
            append!(MK_list,MK_current_best)
            println("MK:",MK_current_best)
        catch
            println("didn't work")
            Number_of_samples = Number_of_samples + 20
            if Number_of_samples > 1000
                break
            end
            continue
        end
        println("p_0_var:",p_0_var)

        if N == 1
            df = DataFrame(Parameter = vcat(first.(p_0_var),Num("Cost")),
                        p_0_value_1 = vcat(last.(p_0_var),MK_current_best)) 
            CSV.write("Intermediate/p_0_values_from_auto.csv",df)
        else
            df = CSV.read("Intermediate/p_0_values_from_auto.csv", DataFrame)
            df[!,Symbol("p_0_value_$N")] = vcat(last.(p_0_var),MK_current_best)
            CSV.write("Intermediate/p_0_values_from_auto.csv",df)
        end

    end
    println("All MK:",MK_list)

end

main()

################################################################################################################################
