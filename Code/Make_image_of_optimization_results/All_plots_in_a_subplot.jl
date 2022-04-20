include("plot_for_all_p.jl")
include("Histogram_plots_of_cost.jl")
include("../../Results/Parameter_values/Results_from_optimization.jl")
plot1 = plot_for_all_p()
plot2 = Histograms_of_cost()

p3 = plot([0 0 0], showaxis = false, grid = false, legend = false)
p2 = plot([0 0 0], showaxis = false, grid = false, legend = false)
p1 = plot([0 0 0], color = [Color_alt_1 Color_alt_2 Color_jalihal], showaxis = false, grid = false, labels = [label_alt_1 label_alt_2 label_jalihal],labelsize=100, legendfontsize=50)
p4 = plot([0 0 0], showaxis = false, grid = false, legend = false)
p5 = plot([0 0 0], showaxis = false, grid = false, legend = false)

plot3 = plot(p1, p2, p3, p4, p5, layout = (5, 1), titlefontsize = 14, guidefontsize = 10, guide_position = :left, margin= 20Plots.mm, foreground_color_text = "white")
plot!(plot3,size=(1500,2500))

leg = false
if leg == false
    plot4 = plot(plot1, plot2, layout = (1, 2), titlefontsize = 14, guidefontsize = 10, guide_position = :left, margin= 20Plots.mm, markersize = 5.5, markercolor="blue", markerstrokewidth=0.2)
    plot!(plot4,size=(3000,2500))
else
    plot4 = plot(plot3, plot1, plot2, layout = (1, 3), titlefontsize = 14, guidefontsize = 10, guide_position = :left, margin= 20Plots.mm, markersize = 5.5, markercolor="blue", markerstrokewidth=0.2)
    plot!(plot4,size=(3750,2500))
end

display(plot4)
savefig(plot4, pwd()*"/Results/Plots_for_new_and_old_paramter_values/plot_and_histogram_for_p_sets.png")

