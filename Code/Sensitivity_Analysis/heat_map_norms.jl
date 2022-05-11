using Plots
using Plots.PlotMeasures
include("sensitivity_analysis.jl")

"""
Plotting heat map of sensitivity vector norms.
"""

abs_senses = hcat(sens_vector_norms_J, sens_vector_norms_1, sens_vector_norms_2)

x_labels = p_var_lookup_table

y_labels = ["Jalihal","Parametervektor 1","Parametervektor 2"]
width = 1000
height = 300


Plots.heatmap(abs_senses',xticks = (1:81,x_labels),yticks = (1:3,y_labels),xrotation=90,size = (width,height),bottom_margin = 75px)