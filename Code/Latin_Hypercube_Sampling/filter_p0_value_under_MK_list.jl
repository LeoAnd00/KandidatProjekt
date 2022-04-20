using CSV
using DataFrames
include("../Model/ODE_functions.jl")
include("../Model/parameter_values.jl") 
function Check_if_p_MK_under_x() 
    df = CSV.read("Intermediate/p_0_values_MK_under_x.csv", DataFrame)

    i = 2 #Choose manually
    df = CSV.read("Intermediate/p_0_values_MK_under_x.csv", DataFrame)
    
    MK_that_works = df[:, Symbol("p_0_value_$i")][1:end]
    #println("MK: ",MK_that_works)

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

Check_if_p_MK_under_x()