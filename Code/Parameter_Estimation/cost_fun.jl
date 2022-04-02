include("../Solve_Jalihal_ODE/Time_course/ODE_methods.jl")
include("../Solve_Jalihal_ODE/Time_course/ODE_functions.jl")
include("../Solve_Jalihal_ODE/Time_course/exp_data.jl")
include("../Solve_Jalihal_ODE/Time_course/parameter_values.jl")

function calc_cost(p_var, p_const)

    function Calc_cost_Glucose_addition_Mig1(p_var, p_const)

        data_for_ode = data__Mig1_glucose_relief
        timelist_for_ode = t_Mig1_glucose_relief
        tspan = (0.0, 20.0) # [min]

        u0_SS = Steady_state_solver(p_const, p_var, (ATP => 0.0, Carbon => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
        
        sol = ODE_solver_FWD(u0_SS, (ATP => 1.0, Carbon => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
        
        u_approx = getindex.(sol.u, 20) #ODE for MIG1 on row 20 in ODE system

        MK_Glucose_addition_Mig1 = sum((data_for_ode - (log.(10,u_approx./(1e-5 .+ 1.0 .- u_approx)))).^2)
    
        return MK_Glucose_addition_Mig1
    end

    function Calc_cost_Glucose_addition_Sch9(p_var, p_const)

        data_for_ode = Minmaxnorm(data_sch9Delta_cAMP, 0.052, 1.227)
        timelist_for_ode = t_sch9Delta_cAMP
        tspan = (0.0, 3.5) # [min]
        # Mutant Sch9_Delta => Sch9_T = 0
        p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 0.0
    
        # Pre-shift => ATP, Carbon = 0
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift => ATP, Carbon = 1
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 6) #ODE for cAMP on row 6 in ODE system

        MK_Glucose_addition_Sch9 = sum((data_for_ode - u_approx).^2)

        p_const[Get_index(p_const_lookup_table, "Sch9_T")] = Sch9_T => 1.0
    
        return MK_Glucose_addition_Sch9
    end

    function Calc_cost_Glucose_addition_cAMP(p_var, p_const)

        data_for_ode = Minmaxnorm(data_cAMP, 0.052, 1.227)
        timelist_for_ode = t_cAMP
        tspan = (0.0, 3.5) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift, Glucose addition => ATP, Carbon = 1
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 6) #ODE for cAMP on row 6 in ODE system
            
        MK_Glucose_addition_cAMP = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_addition_cAMP
    end

    function Calc_cost_Glucose_addition_Sch9_p(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Sch9_glucose_relief)
        timelist_for_ode = t_Sch9_glucose_relief
        tspan = (0.0, 30) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system

        MK_Glucose_addition_Sch9_p = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_addition_Sch9_p
    end

    function Calc_cost_Glucose_starvation_Snf1_p(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Snf1)
        timelist_for_ode = t_Snf1
        tspan = (0.0, 62) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift, Glucose starvation => Carbon, ATP = 0
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
        u_approx = getindex.(sol.u, 10) #ODE for Snf1 on row 10 in ODE system

        MK_Glucose_starvation_Snf1_p = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_starvation_Snf1_p
    end

    function Calc_cost_Glucose_starvation_Sch9_p(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Sch9_glucose_starve)
        timelist_for_ode = t_Sch9_glucose_starve
        tspan = (0.0, 30) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system

        MK_Glucose_starvation_Sch9_p = sum((data_for_ode - u_approx).^2)

        return MK_Glucose_starvation_Sch9_p
    end

    function Calc_cost_Sch9P_glutamine_L(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Sch9P_glutamine_L, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_L
        tspan = (0.0, 32) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady state för parametrarna p
    
        # Low glutamine => Glutamine_ext = 0.3
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system

        MK_Sch9P_glutamine_L = sum((data_for_ode - u_approx).^2)

        return MK_Sch9P_glutamine_L
    end

    function Calc_cost_Sch9P_glutamine_H(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Sch9P_glutamine_H, 4.82, 52.57)
        timelist_for_ode = t_Sch9P_glutamine_H
        tspan = (0.0, 32) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0)) # Returnerar steady state för parametrarna p
    
        # High Glutamine => Glutamine_ext = 1.0
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system
 
        MK_Sch9P_glutamine_H = sum((data_for_ode - u_approx).^2)

        return MK_Sch9P_glutamine_H
    end

    function Calc_cost_Glutamine_addition_Sch9_gtr1Delta(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Sch9_gtr1Delta, 4.82, 52.27)
        timelist_for_ode = t_Sch9_gtr1Delta
        tspan = (0.0, 32) # [min]
    
        w_torc_ego_true = p_var[Get_index(p_var_lookup_table, "w_torc_ego")]
        w_torc_egoin_true = p_var[Get_index(p_var_lookup_table, "w_torc_egoin")]
    
        # Mutant gtr1_Delta => EGO_T = 0, w_torc_ego = 0, w_torc_egoin = 0 i 
        p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 0.0
        p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego => 0.0
        p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin => 0.0
    
        # Pre-shift => Glutamine_ext = 0, ATP, Carbon = 1 
        u0_SS = Steady_state_solver(p_const, p_var, (Glutamine_ext => 0.0, Carbon => 1.0, ATP => 1.0)) 
    
        # High glutamine => Glutamine_ext = 1
        sol = ODE_solver_FWD(u0_SS, (Glutamine_ext => 1.0, Carbon => 1.0, ATP => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 12) #ODE for Sch9 on row 12 in ODE system
                
        MK_Glutamine_addition_Sch9_gtr1Delta = sum((data_for_ode - u_approx).^2)
    
        p_const[Get_index(p_const_lookup_table, "EGO_T")] = EGO_T => 1.0
        p_var[Get_index(p_var_lookup_table, "w_torc_ego")] = w_torc_ego_true
        p_var[Get_index(p_var_lookup_table, "w_torc_egoin")] = w_torc_egoin_true
    
        return MK_Glutamine_addition_Sch9_gtr1Delta
    end

    function Calc_cost_Rapamycin_treatment(p_var, p_const)

        data_for_ode = Minmaxnorm(data_Rib_rap)
        timelist_for_ode = t_Rib_rap
        tspan = (0.0, 92.0) # [min]
    
        u0_SS = Steady_state_solver(p_const, p_var, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)) # Returnerar steady state för parametrarna p
    
        # Post-shift: Rapamycin treatment => TORC1_T = 0.0 
        p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 0.0
    
        sol = ODE_solver_FWD(u0_SS, (Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0), tspan, p_const, p_var, timelist_for_ode)
    
        u_approx = getindex.(sol.u, 25) #ODE for RIB on row 25 in ODE system

        MK_Rapamycin_treatment = sum((data_for_ode - (u_approx./(1e-3 .+ u0_SS[25]))).^2)
    
        p_const[Get_index(p_const_lookup_table, "TORC1_T")] = TORC1_T => 1.0
    
        return MK_Rapamycin_treatment
    end

    MK_Glucose_addition_Mig1 = Calc_cost_Glucose_addition_Mig1(p_var, p_const)
    MK_Glucose_addition_Sch9 = Calc_cost_Glucose_addition_Sch9(p_var, p_const)
    MK_Glucose_addition_cAMP = Calc_cost_Glucose_addition_cAMP(p_var, p_const)
    MK_Glucose_addition_Sch9_p = Calc_cost_Glucose_addition_Sch9_p(p_var, p_const)
    MK_Glucose_starvation_Snf1_p = Calc_cost_Glucose_starvation_Snf1_p(p_var, p_const)
    MK_Glucose_starvation_Sch9_p = Calc_cost_Glucose_starvation_Sch9_p(p_var, p_const)
    MK_Sch9P_glutamine_L = Calc_cost_Sch9P_glutamine_L(p_var, p_const)
    MK_Sch9P_glutamine_H = Calc_cost_Sch9P_glutamine_H(p_var, p_const)
    MK_Glutamine_addition_Sch9_gtr1Delta = Calc_cost_Glutamine_addition_Sch9_gtr1Delta(p_var, p_const)
    MK_Rapamycin_treatment = Calc_cost_Rapamycin_treatment(p_var, p_const)

    MK = MK_Glucose_addition_Mig1 + MK_Glucose_addition_Sch9 + MK_Glucose_addition_cAMP + MK_Glucose_addition_Sch9_p + MK_Glucose_starvation_Snf1_p + MK_Glucose_starvation_Sch9_p + MK_Sch9P_glutamine_L + MK_Sch9P_glutamine_H + MK_Glutamine_addition_Sch9_gtr1Delta + MK_Rapamycin_treatment

    #prob = ODEProblem(ODE_sys,init_vals,tspan,params)
    #sol = solve(prob,Rodas5(),saveat=tvals,verbose = false)
    #if sol.retcode != :Success
    #    return Inf
    #end
    #return cost_fun(sol.u)
    return MK
end