"""
    plot_for_all_ODE_solutions_with_data()

    Makes a subplot for all ODE solutions with data for all 3 parametervectors.
"""
function plot_for_all_ODE_solutions_with_data()
    include("Glucose_addition_Mig1.jl")
    include("Glucose_addition_Sch9.jl")
    include("Glucose_addition.jl")
    include("Glucose_starvation.jl")
    include("Glutamine_addition_gtr1.jl")
    include("Glutamine_addition.jl")
    include("Rapamycin_treatment.jl")

    plot_Glucose_addition_Mig1 = p1 = Glucose_addition_Mig1()
    plot_Glucose_addition_Sch9_cAMP = p2 = Glucose_addition_Sch9_cAMP()
    plot_Glucose_addition_cAMP = p3 = Glucose_addition_cAMP()
    plot_Glucose_addition_Sch9 = p4 = Glucose_addition_Sch9()
    plot_Glucose_starvation_Snf1 = p5 = Glucose_starvation_Snf1()
    plot_Glucose_starvation_Sch9 = p6 = Glucose_starvation_Sch9()
    plot_Glutamine_addition_gtr1 = p7 = Glutamine_addition_gtr1()
    plot_Glutamine_addition_L = p8 = Glutamine_addition_L()
    plot_Glutamine_addition_H = p9 = Glutamine_addition_H()
    plot_Rapamycin_treatment = p10 = Rapamycin_treatment()

    """
    display(plot_Glucose_addition_Mig1)
    display(plot_Glucose_addition_Sch9_cAMP)
    display(plot_Glucose_addition_cAMP)
    display(plot_Glucose_addition_Sch9)
    display(plot_Glucose_starvation_Snf1)
    display(plot_Glucose_starvation_Sch9)
    display(plot_Glutamine_addition_gtr1)
    display(plot_Glutamine_addition_L)
    display(plot_Glutamine_addition_H)
    display(plot_Rapamycin_treatment)
    """

    plot1 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, layout = (5, 2), legend = false, titlefontsize = 14, guidefontsize = 10, guide_position = :left, margin= 20Plots.mm)
    plot!(plot1,size=(1500,2500))
    #display(plot1)
    return plot1
end