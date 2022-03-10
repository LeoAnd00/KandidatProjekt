include("ODE_methods.jl")

function Steady_state_pertubation_solver(model_inputs_preshift, model_inputs_postshift, output_variable)
    include("parameter_values.jl")

    # Pre-shift
    u0_SS_preshift = Steady_state_solver(p, model_inputs_preshift) # Returnerar steady state för parametrarna p

    # Post-shift
    u0_SS_postshift = Steady_state_solver(p, model_inputs_postshift)

    index = findfirst(x->x==string(output_variable), u_lookup_table)

    println("Pre-shift: ", u0_SS_preshift[index], ", post-shift: ", u0_SS_postshift[index])
end

function Steady_state_pertubation_solver(mutations, model_inputs_preshift, model_inputs_postshift, output_variable)
    include("parameter_values.jl")

    # Implement mutations
    for i in range(1, length(mutations))
        p[findfirst(x->x==string(first.(mutations)[i]), p_lookup_table)] = last.(mutations)[i]
    end

    # Pre-shift
    u0_SS_preshift = Steady_state_solver(p, model_inputs_preshift) # Returnerar steady state för parametrarna p

    # Post-shift
    u0_SS_postshift = Steady_state_solver(p, model_inputs_postshift)

    index = findfirst(x->x==string(output_variable), u_lookup_table)

    println("Pre-shift: ", u0_SS_preshift[index], ", post-shift: ", u0_SS_postshift[index])
end