using CSV
using DataFrames
using LinearAlgebra
using DifferentialEquations
using Plots
using DiffEqSensitivity
using SteadyStateDiffEq
using ForwardDiff
using ProgressMeter
using ModelingToolkit
include("../Model/ODE_functions.jl") 
include("../Model/parameter_values.jl")
include("../Model/ODE_methods.jl")
include("../../Data/exp_data.jl")
include("Optimization.jl")

"""
    main_Op_p_0()

    Optimize the cost for multiple p_0_var values and puts the result in a CSV file
"""
function main_Op_p_0()
    include("../Model/parameter_values.jl")
    
    df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
    List_of_p_that_didnt_work = []

    for i in 1:(length(df[1,:])-1)
        println(" I : ",i)
        
        p_names_var = first.(p_var)

        df = CSV.read("Intermediate/p_0_values.csv", DataFrame)

        p_0_var = df[:, Symbol("p_0_value_$i")][1:end-1]

        p_var_temp = [Pair(p_names_var[1], p_0_var'[1])]
        for item in range(2, length(p_0_var))
            B = Pair(p_names_var[item], p_0_var'[item])
            p_var_temp = vcat(p_var_temp,B)
        end
        p_0_var = p_var_temp

        MK = []
        p_var_result = []
        try
            MK, p_var_result = Optimize_p(p_0_var)
        catch
            append!(List_of_p_that_didnt_work,i)
            continue
        end

        parameters_and_MK = vcat(first.(p_var),Num("Cost"))

        values_of_p_and_MK = vcat(last.(p_var_result),MK[end])
        Prior_values_of_p_and_MK = vcat(last.(p_0_var),MK[1])
        #print("length: ",length(Prior_values_of_p_and_MK))
        if i == 1
            df = DataFrame(Parameter = parameters_and_MK,
                        p_0_value_1 = Prior_values_of_p_and_MK,
                        Optimum_from_p_0_value_1 = values_of_p_and_MK) 
            CSV.write("Intermediate/optimization_of_p0.csv",df)
        else
            df = CSV.read("Intermediate/optimization_of_p0.csv", DataFrame)
            #print(" length_p_0_v: ",length(df[:, Symbol("Optimum_from_p_0_value_1")]))
            df[!,Symbol("p_0_value_$i")] = Prior_values_of_p_and_MK
            df[!,Symbol("Optimum_from_p_0_value_$i")] = values_of_p_and_MK
            CSV.write("Intermediate/optimization_of_p0.csv",df)
        end
    end

    println("List_of_p_that_didnt_work: ",List_of_p_that_didnt_work)
end

"""
    Optimize_p(p_0_var)

    Optimize the cost for p_0_var.
"""
function Optimize_p(p_0_var)
    MK, p_var_reslut = main_FWD_optimization(p_0_var)
    return MK, p_var_reslut
end

main_Op_p_0()