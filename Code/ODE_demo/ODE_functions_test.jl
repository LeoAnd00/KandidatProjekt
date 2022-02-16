#Definerar ODE- systemes 3 ekvationer
# u[1]=x, u[2]=y, u[3]=z
include("Get_parmeter_values_test.jl")

function ODE_test(du, u,p, t)
    du[1] = dict1["P1"]*(u[2]-u[1])
    du[2] = u[1]*(dict1["P2"]-u[3]) - u[2]
    du[3] = u[1]*u[2] - dict1["P3"]*u[3]
end