using CSV
using DataFrames
include("../Model/ODE_functions.jl")
include("../Model/parameter_values.jl") 

"""
    Get_one_vector_from_p0_values_list(x) 

    Gives the user the parametervector in a form suitable for ModelingToolkit 
    for parametervector i in p_0_values.csv that in the Intermediate folder.
"""
function Get_one_vector_from_p0_values_list()
    df = CSV.read("Intermediate/p_0_values.csv", DataFrame)
    i = 1 # The vector that will be displayed

    p_names_var = first.(p_var)

    p_0_var = df[:, Symbol("p_0_value_$i")][1:end-1]

    p_var_temp = [Pair(p_names_var[1], p_0_var'[1])]
    for item in range(2, length(p_0_var))
        B = Pair(p_names_var[item], p_0_var'[item])
        p_var_temp = vcat(p_var_temp,B)
    end
    p_0_var = p_var_temp
    println("p 0 var:",p_0_var)
    println("MK: ",df[:, Symbol("p_0_value_$i")][end])
end

Get_one_vector_from_p0_values_list()