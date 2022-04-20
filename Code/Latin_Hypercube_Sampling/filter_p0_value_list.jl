using CSV
using DataFrames
include("../Model/ODE_functions.jl")
include("../Model/parameter_values.jl") 
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

x = 1.6
Check_if_p_MK_under_x(x)