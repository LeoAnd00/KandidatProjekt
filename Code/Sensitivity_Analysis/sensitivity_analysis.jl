using LinearAlgebra
include("../../Data/exp_data.jl")
include("../Model/ODE_methods.jl")
include("../Model/parameter_values.jl")
include("../Model/ODE_functions.jl")
include("shift_info_SA.jl")
include("../../Results/Parameter_values/Results_from_optimization.jl")

function dy_dp(param,out_index,i)
    sim_param = param
    pre_p_const = p_const
    for (p_key,value) in all_pre_p_const_changes[i]
        pre_p_const[Get_index(p_const_lookup_table,p_key)] = value
    end
    post_p_const = p_const
    for (p_key,value) in all_post_p_const_changes[i]
        post_p_const[Get_index(p_const_lookup_table,p_key)] = value
    end
    for (p_key,value) in all_pre_p_var_changes[i]
        sim_param[Get_index(p_var_lookup_table,p_key)] = value
    end

    u0 = Steady_state_solver(last.(pre_p_const),sim_param,pre_shifts[i])

    sens_sol = sensitivity_solver(u0,post_shifts[i],all_tvals[i],last.(post_p_const),sim_param)    
    if sens_sol.retcode != :Success
        throw(ErrorException("instability"))
    end

    sens_mat = extract_local_sensitivities(sens_sol)

    dy_dp = []
    for t in 1:length(all_tvals[i])
        push!(dy_dp, [])
        for p in 1:length(param)
            push!(dy_dp[t], sens_mat[2][p][out_index,t])
        end
    end    

    return dy_dp
end

function sens_mat(p)
    sens_mat = []
    for i in 1:length(out_indices)
        block = dy_dp(p,out_indices[i],i)
        append!(sens_mat,block)
    end
    return sens_mat
end

function cond_nr(sens)
    sm = zeros(81,86)
    for i in 1:86
        sm[:,i] = sens[i]
    end
    return cond(sm)
end

sens_Jalihal = sens_mat(last.(p_var))
sens_alt1 = sens_mat(last.(p_var_alt_1))
sens_alt2 = sens_mat(last.(p_var_alt_2))
sens_alt3 = sens_mat(last.(p_var_alt_3))

cond_nr_J = cond_nr(sens_Jalihal)
cond_nr_1 = cond_nr(sens_alt1)
cond_nr_2 = cond_nr(sens_alt2)
cond_nr_3 = cond_nr(sens_alt3)