using CSV
using DataFrames
include("../Model/ODE_functions.jl")
include("../Model/parameter_values.jl") 

"""
    filter_p0_values_opt_list() 

    Gives the user the parametervector in a form suitable for ModelingToolkit 
    for parametervector i in optimization_of_p0.csv.
"""
function filter_p0_values_opt_list()
    df = CSV.read("Intermediate/optimization_of_p0.csv", DataFrame)

    i = 11 #Choose manually

    p_names_var = first.(p_var)

    p_0_var = df[:, Symbol("Optimum_from_p_0_value_$i")][1:end-1]

    p_var_temp = [Pair(p_names_var[1], p_0_var'[1])]
    for item in range(2, length(p_0_var))
        B = Pair(p_names_var[item], p_0_var'[item])
        p_var_temp = vcat(p_var_temp,B)
    end
    p_0_var = p_var_temp
    println("p 0 var:",p_0_var)
    println("MK: ",df[:, Symbol("Optimum_from_p_0_value_$i")][end])
end

filter_p0_values_opt_list()