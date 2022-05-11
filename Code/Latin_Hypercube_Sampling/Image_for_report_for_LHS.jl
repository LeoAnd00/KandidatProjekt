using LatinHypercubeSampling
using Plots

# Creates a figure that illustrates the difference between LHS and random 
# sampling

plan = randomLHC(20,2)
scaled_plan = scaleLHC(plan,[(0.0,1.0),(0.0,1.0)]) # Defines the interval 
# for which values the 2 dimensions can take

gr() 
plotd = plot(scaled_plan[:,1], scaled_plan[:,2], seriestype = :scatter, 
title = "Latin Hypercube Sampling", legend = false)

plotd2 = plot(rand(20), rand(20), seriestype = :scatter, 
title = "Random Sampling", legend = false, color = "orange")
l = @layout [a b]

plot_result = plot(plotd, plotd2,layout = l, gridalpha=0.3,xticks=0.0:0.2:1.0,
yticks=0.0:0.2:1.0,tickfontcolor = "black", minorgrid = true, minorticks = 4, 
minorgridalpha = 0.3, aspect_ratio=:equal, markersize = 3)
xlabel!("Dimension 1")
ylabel!("Dimension 2")
xlims!((-0.1, 1.1))
ylims!((-0.1, 1.1))

display(plot_result)

Want_to_plot = false
if Want_to_plot == true
    savefig(plot_result, pwd()*"/Results/LHS/LHS_png.png")
    savefig(plot_result, pwd()*"/Results/LHS/LHS_svg.svg")
end
