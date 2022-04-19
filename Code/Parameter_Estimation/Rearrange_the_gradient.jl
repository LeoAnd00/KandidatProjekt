
"""
    Rearrange_the_gradient(p_var_FWD, p_const_FWD, dMK_dp)

    Rearrange the gradient so that they are in the correct order.
"""
function Rearrange_the_gradient(p_var_FWD, p_const_FWD, dMK_dp)

    tspan = (0.0, 30.0) # [min]
    u0_SS = ones(1,28) 

    p_var_FWD_new = []
    for i in 1:length(last.(p_var_FWD))
        append!(p_var_FWD_new,1 + 0.01*i)
    end

    p_const_FWD_new = []
    for i in 1:length(last.(p_const_FWD))
        append!(p_const_FWD_new,2 + 0.01*i)
    end

    p_values_process = vcat(p_var_FWD_new, p_const_FWD_new)

    p_original = vcat(p_var_FWD,p_const_FWD)

    p_values_dual_temp = [Pair(first.(p_original)[1], p_values_process'[1])]
    for i in range(2, length(last.(p_original)))
        B = Pair(first.(p_original)[i], p_values_process'[i])
        p_values_dual_temp = vcat(p_values_dual_temp,B)
    end
    p_new = p_values_dual_temp

    prob = ODEProblem(ODE_sys, u0_SS, tspan, p_new)

    all_p_values_after_prob = prob.p

    p_var_FWD_values_after_prob = []
    for i in all_p_values_after_prob
        if i < 2.0
            append!(p_var_FWD_values_after_prob,round((i-1)*100; digits = 3))
        end
    end

    Pre_process_grad = []

    for i in 1:length(dMK_dp)
        if all_p_values_after_prob[i] < 2
            append!(Pre_process_grad,dMK_dp[i])
        end
    end

    finished_process_grad = []
    for i in 1:length(p_var_FWD_values_after_prob)
        N = 0
        for item in p_var_FWD_values_after_prob
            N += 1
            if item == i
                append!(finished_process_grad,Pre_process_grad[N])
                break
            end
        end
    end

    return finished_process_grad
end




