using Plots
using Plots.PlotMeasures
include("sensitivity_analysis.jl")

"""
Plotting heat map of pairwise sensitivity vector dependencies.
"""

x_labels = p_var_lookup_table
y_labels = x_labels
width = 1000
height = 1000

Plots.heatmap(rel_J,xticks = (1:81,x_labels),yticks = (1:81,y_labels),xrotation=90,
size = (width,height),bottom_margin = 75px)
Plots.heatmap(rel_1,xticks = (1:81,x_labels),yticks = (1:81,y_labels),xrotation=90,
size = (width,height),bottom_margin = 75px)
Plots.heatmap(rel_2,xticks = (1:81,x_labels),yticks = (1:81,y_labels),xrotation=90,
size = (width,height),bottom_margin = 75px)