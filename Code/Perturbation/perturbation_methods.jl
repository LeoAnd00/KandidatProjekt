include("../Model/ODE_methods.jl")

function Steady_state_pertubation_solver(model_inputs_preshift, model_inputs_postshift, output_variable)
    include("../Model/parameter_values.jl")

    # Pre-shift
    u0_SS_preshift = Steady_state_solver(p_const, p_var, model_inputs_preshift) # Returnerar steady state för parametrarna p

    # Post-shift
    u0_SS_postshift = Steady_state_solver(p_const, p_var, model_inputs_postshift)

    index = Get_index(u_lookup_table, string(output_variable))

    return (u0_SS_preshift[index] => u0_SS_postshift[index])
end

function Steady_state_pertubation_solver(p_changes, model_inputs_preshift, model_inputs_postshift, output_variable)
    include("../Model/parameter_values.jl")

    # Implement parameter changes
    for i in range(1, length(p_changes))
        p_conc[Get_index(p_conc_lookup_table, string(first.(p_changes)[i]))] = first.(p_changes)[i] => last.(p_changes)[i]
    end

    # Pre-shift
    u0_SS_preshift = Steady_state_solver(p_conc, model_inputs_preshift) # Returnerar steady state för parametrarna p

    # Post-shift
    u0_SS_postshift = Steady_state_solver(p_conc, model_inputs_postshift)

    index = Get_index(u_lookup_table, string(output_variable))

    return (u0_SS_preshift[index] => u0_SS_postshift[index])
end

function Minmaxnorm_shifts(shift)

    for i in 1:length(shift)
        c = shift[i]
        min = minimum([c.first, c.second])
        max = maximum([c.first, c.second])

        norm_vals = Minmaxnorm([c.first, c.second], min, max)
        shift[i] = norm_vals[1] => norm_vals[2]
    end

    return shift
end

function Minmaxnorm_shifts(shift, min, max)

    for i in 1:length(shift)
        c = shift[i]
        norm_vals = Minmaxnorm([first.(c) last.(c)], min, max)
        shift[i] = norm_vals[1] => norm_vals[2]
    end

    return shift
end

function Minmaxnorm_data(data)

    min = minimum([first.(data) last.(data)])
    max = maximum([first.(data) last.(data)])
    
    for i in 1:length(data)
        c = data[i]
        norm_vals = Minmaxnorm([first.(c) last.(c)], min, max)
        data[i] = norm_vals[1] => norm_vals[2]
    end

    return data
end

function Minmaxnorm_data(data, min, max)
    
    for i in 1:length(data)
        c = data[i]
        norm_vals = Minmaxnorm([first.(c) last.(c)], min, max)
        data[i] = norm_vals[1] => norm_vals[2]
    end

    return data
end