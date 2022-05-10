using LinearAlgebra
include("../../Data/exp_data.jl")
include("../Model/ODE_methods.jl")
include("../Model/parameter_values.jl")
include("../Model/ODE_functions.jl")
include("shift_info_SA.jl")
include("../../Results/Parameter_values/Results_from_optimization.jl")

"""
    dy_dp(param,out_index,i)

    Compute the sensitivities of out signal out_index in the point param.
"""
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

"""
    sens_mat(p)

    Compute the sensitivity matrix in the point p.
"""
function sens_mat(p)
    sens_mat = []
    for i in 1:length(out_indices)
        block = dy_dp(p,out_indices[i],i)
        append!(sens_mat,block)
    end
    sm = zeros(81,86)
    for i in 1:86
        sm[:,i] = sens_mat[i]
    end
    return sm
end

"""
    sum_sens_vecs(sens_mat)

    Compute the sensitivity vector norms and the sum of them for sensitivity matrix sens_mat.
"""
function sum_vec_sens(sens_mat)
    s = 0
    v = []
    for i in 1:length(sens_mat[:,1])
        n = norm(sens_mat[i,:])
        s += n
        append!(v,n)
    end
    return s,v
end


"""
    cos_rel_senses(sens_mat)

    Compute a matrix of cosines of the angles between all pairs of
    sensitivity vectors in sensitivity matrix sens_mat.
"""
function cos_rel_senses(sens_mat)
    m = zeros(81,81)
    for i in 1:81
        for j in 1:81
            n = norm(sens_mat[i,:])*norm(sens_mat[j,:])
            if n == 0.0
                m[i,j] = -1.0
            else
                m[i,j] = abs(dot(sens_mat[i,:],sens_mat[j,:])/n)
            end
        end
    end
    return m
end

sens_Jalihal = sens_mat(last.(p_var))
sens_alt1 = sens_mat(last.(p_var_alt_1))
sens_alt2 = sens_mat(last.(p_var_alt_2))

cond_nr_J = cond(sens_Jalihal)
cond_nr_1 = cond(sens_alt1)
cond_nr_2 = cond(sens_alt2)

rank_J = rank(sens_Jalihal)
rank_1 = rank(sens_alt1)
rank_2 = rank(sens_alt2)

abs_sens_J,sens_vector_norms_J = sum_vec_sens(sens_Jalihal)
abs_sens_1,sens_vector_norms_1 = sum_vec_sens(sens_alt1)
abs_sens_2,sens_vector_norms_2 = sum_vec_sens(sens_alt2)

rel_J = cos_rel_senses(sens_Jalihal)
rel_1 = cos_rel_senses(sens_alt1)
rel_2 = cos_rel_senses(sens_alt2)

out = "k_J = "*string(cond_nr_J)*"\n
k_1 = "*string(cond_nr_1)*"\n
k_2 = "*string(cond_nr_2)*"\n
r_J = "*string(rank_J)*"\n
r_1 = "*string(rank_1)*"\n
r_2 = "*string(rank_2)*"\n
sum_sens_vecs J = "*string(abs_sens_J)*"\n
sum_sens_vecs 1 = "*string(abs_sens_1)*"\n
sum_sens_vecs 2 = "*string(abs_sens_2)*"\n
sens_vecs_norms J = "*string(sens_vector_norms_J)*"\n
sens_vecs_norms 1 = "*string(sens_vector_norms_1)*"\n
sens_vecs_norms 2 = "*string(sens_vector_norms_2)*"\n
sens_mat_J = "*string(sens_Jalihal)*"\n
sens_mat_1 = "*string(sens_alt1)*"\n
sens_mat_2 = "*string(sens_alt2)*"\n
rel_sens_J = "*string(rel_J)*"\n
rel_sens_1 = "*string(rel_1)*"\n
rel_sens_2 = "*string(rel_2)*"\n
rel_sens_J = "*string(rel_J)*"\n
rel_sens_1 = "*string(rel_1)*"\n
rel_sens_2 = "*string(rel_2);

open("Results/Parameter_values/sensitivities_out.txt","w") do file
    write(file,out)
end