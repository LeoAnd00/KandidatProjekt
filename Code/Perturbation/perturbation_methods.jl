include("../Model/ODE_methods.jl")
include("../Model/parameter_values.jl")

"""
    Steady_state_pertubation_solver(model_inputs_preshift, model_inputs_postshift, output_variable)
    Steady_state_pertubation_solver(p_changes, model_inputs_preshift, model_inputs_postshift, output_variable)

Calculate the steady state values both for pre- and postshift nutrient conditions. Values are returned as a Pair{Float64, Float64}
"""
function Steady_state_pertubation_solver(model_inputs_preshift, model_inputs_postshift, output_variable)

    u0_SS_preshift = Steady_state_solver(p_const, p_var, model_inputs_preshift) 
    u0_SS_postshift = Steady_state_solver(p_const, p_var, model_inputs_postshift)

    index = Get_index(u_lookup_table, string(output_variable))

    return (u0_SS_preshift[index] => u0_SS_postshift[index])
end

function Steady_state_pertubation_solver(p_changes, model_inputs_preshift, model_inputs_postshift, output_variable)

    p_conc_copy = copy(p_conc)

    # implement parameter changes
    for i in range(1, length(p_changes))
        p_conc_copy[Get_index(p_conc_lookup_table, string(first.(p_changes)[i]))] = 
        first.(p_changes)[i] => last.(p_changes)[i]
    end

    u0_SS_preshift = Steady_state_solver(p_conc_copy, model_inputs_preshift) 
    u0_SS_postshift = Steady_state_solver(p_conc_copy, model_inputs_postshift)

    index = Get_index(u_lookup_table, string(output_variable))

    return (u0_SS_preshift[index] => u0_SS_postshift[index])
end

"""
    Minmaxnorm_shifts(shifts_array)
    Minmaxnorm_shifts(shifts_array, min, max)

Normalizes the results of simulations of a type of shift according to `Minmaxnorm()`, using 
either min/max retrived from the shift value or as function inputs.
"""
function Minmaxnorm_shifts(shifts_array)

    for i in 1:length(shifts_array)
        c = shifts_array[i]
        min = minimum([c.first, c.second])
        max = maximum([c.first, c.second])

        norm_vals = Minmaxnorm([c.first, c.second], min, max)
        shifts_array[i] = norm_vals[1] => norm_vals[2]
    end

    return shifts_array
end

function Minmaxnorm_shifts(shift, min, max)

    for i in 1:length(shift)
        c = shift[i]
        norm_vals = Minmaxnorm([first.(c) last.(c)], min, max)
        shift[i] = norm_vals[1] => norm_vals[2]
    end

    return shift
end

"""
    Minmaxnorm_shifts(shifts_array)
    Minmaxnorm_shifts(shifts_array, min, max)

Normalizes the data according to `Minmaxnorm()`, using either min, max retrived from the data or as function inputs 
"""
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