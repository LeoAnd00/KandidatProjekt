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

"""
    Check_if_p_works(p_var, p_const, scaled_plan)

    Checks if the parameter sets from the latin hypercube works and saves them
    into a CSV file.
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

    return N
end

"""
    main()

    First one must define a list called p_0_values_list.csv with 
    prametervectors that the sampling will start around. This list can for 
    example be created by first generating a more general latin hypercube with
    latin_hypercube_sampling.jl and one can then filter the results to get a 
    list with parametervectors that gave a cost below a choosen value. Then 
    the code (main()) will create a latin hypercube based on the 
    Number_of_samples and what intervals the user defines. If it finds a 
    parametervector with a cost below the previous best one it will start 
    generating a new latin hypercube around the new value. If none of the 
    costs were below the previously best one it will try to increase the 
    number of samples by a choosen amount and generate a new latin hypercube. 
    If the cost is less then the previous best but not good enough, specified 
    by how much the user wants the cost to decrease, it will move on to the 
    next parametervector in p_0_values_list.csv. The code will also save all 
    intermediate results in csv files.

    Note that one must do some modifications in 
    Check_if_p_works(p_var, p_const, scaled_plan) if the user wants to give 
    the intervals on logarithmic-scale.
"""
function main()
    include("../Model/parameter_values.jl")
    N_3 = 0
    while true
        N_3 += 1
        MK_current_best = 1e100
        MK_list = []
        p_0_var = p_var
        N = 0
        Number_of_samples = 300
        N_2 = 1
        MK_prev = 0
        while true
            N += 1
            p_var = p_0_var
            if N > 1
                plan = randomLHC(Number_of_samples,81) #samples, dimensions

                plan_scale = [(0.8*last.(p_var)[1],1.2*last.(p_var)[1])]
                for i in 1:80
                    push!(plan_scale,(0.8*last.(p_var)[i+1],
                    1.2*last.(p_var)[i+1]))
                end

                scaled_plan = scaleLHC(plan,plan_scale)

                N_check = Check_if_p_works(p_var, p_const, scaled_plan)
                println(" Number of working parameter sets: ", N_check) 
            end

            if N > 1
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
                    Number_of_samples = Number_of_samples + 200 # Adds this 
                    # amount to number of samples if no cost under the 
                    # previous best was found
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
                    continue
                end
                println("p_0_var:",p_0_var)
            else
                Best_i = N_3

                df = CSV.read("Intermediate/p_0_values_list.csv", DataFrame)
                p_names_var = first.(p_var)
                p_0_var = df[:, Symbol("p_0_value_$Best_i")][1:end-1]
                MK_df = df[:, Symbol("p_0_value_$Best_i")][end]
    
                p_var_temp = [Pair(p_names_var[1], p_0_var'[1])]
                for item in range(2, length(p_0_var))
                    B = Pair(p_names_var[item], p_0_var'[item])
                    p_var_temp = vcat(p_var_temp,B)
                end
                p_0_var = p_var_temp
                MK_current_best = MK_df
                append!(MK_list,MK_current_best)
                println("MK:",MK_current_best)

                println("p_0_var:",p_0_var)
            end

            if N == 1
                df = DataFrame(Parameter = vcat(first.(p_0_var),Num("Cost")),
                            p_0_value_1 = vcat(last.(p_0_var),MK_current_best)) 
                CSV.write("Intermediate/p_0_values_from_auto_$N_3.csv",df)
            else
                N_2 += 1
                df = CSV.read("Intermediate/p_0_values_from_auto_$N_3.csv"
                , DataFrame)
                df[!,Symbol("p_0_value_$N_2")] = vcat(last.(p_0_var),
                MK_current_best)
                CSV.write("Intermediate/p_0_values_from_auto_$N_3.csv",df)
            end

            if N <= 1
                MK_prev = 1e10 # Just to make sure grad_of_cost isn't less 
                # then 0.05 for the first iteration
            end

            grad_of_cost = abs(MK_df - MK_prev)
            if grad_of_cost < 0.05 # The value of which the cost must 
                # decrease from each iteration
                break
            end

            if N > 1
                MK_prev = MK_df # Used in next iteration
            end

        end
    end
    println("All MK:",MK_list)

end

main()
