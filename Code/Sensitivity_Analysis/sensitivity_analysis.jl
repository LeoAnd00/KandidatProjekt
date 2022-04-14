include("../Solve_Jalihal_ODE/Time_course/exp_data.jl")
include("../Solve_Jalihal_ODE/Time_course/ODE_methods.jl")
include("../Solve_Jalihal_ODE/Time_course/parameter_values.jl")
include("../Solve_Jalihal_ODE/Time_course/ODE_functions.jl")
include("shifts.jl")

function sim_out(param,key)
    sim_param = param

    tvals = get(tvals_dict,key,nothing)
    if isnothing(tvals)
        throw(KeyError(key))
    end
    pre_shift = get(pre_shift_dict,key,nothing)
    post_shift = get(post_shift_dict,key,nothing)
    pre_p_const = p_const
    for (p_key,value) in get(pre_p_const_changes,key,[])
        pre_p_const[Get_index(p_const_lookup_table,p_key)] = value
    end
    post_p_const = p_const
    for (p_key,value) in get(post_p_const_changes,key,[])
        post_p_const[Get_index(p_const_lookup_table,p_key)] = value
    end
    for (p_key,value) in get(pre_p_var_changes,key,[])

        sim_param[Get_index(p_var_lookup_table,p_key)] = value
    end

    
    u0 = Steady_state_solver(last.(pre_p_const),sim_param,pre_shift)

    sens_sol = sensitivity_solver(u0,post_shift,tvals,last.(post_p_const),sim_param)    
    if sens_sol.retcode != :Success
        throw(ErrorException("instability"))
    end

    sens_mat = extract_local_sensitivities(sens_sol)

    u_key = get(u_dict,key,key)
    u = Get_index(u_lookup_table,u_key*"(t)")
    # y = []
    dy_dp = []

    for t in 1:length(tvals)
        # push!(y, sens_sol.u[t][u])
        push!(dy_dp, [])
        for p in 1:length(param)
            push!(dy_dp[t], sens_mat[2][p][u,t])
        end
    end    

    return dy_dp
end

function get_sensitivity_matrix(p)
    sensitivity_matrix = Dict()
    for (key,_) in dat_dict
        try
            dy_dp = sim_out(p,key)
            sensitivity_matrix[key] = dy_dp
        catch e
            throw(ErrorException(e.msg))
        end
    end
    return sensitivity_matrix
end