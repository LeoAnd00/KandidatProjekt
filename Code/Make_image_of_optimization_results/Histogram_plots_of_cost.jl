using StatsPlots
function Histograms_of_cost()
  include("../../Results/Parameter_values/Results_from_optimization.jl")
  #include("Cost_for_histograms.jl")
  #Cost_alt_1, Cost_alt_2, Cost_jalihal = Cost_for_histograms()
  #println(Cost_alt_1)
  #println(Cost_alt_2)
  #println(Cost_jalihal)
  Cost_alt_1 = [0.02771829851138364, 0.11530809012155227, 0.14658811799803514, 0.038428681050842434, 0.02894780772394849, 0.1127683119230797, 0.06509447763039686, 0.10031530561560253, 0.08076111254137344, 0.0036672576122928754]
  Cost_alt_2= [0.05048587297214781, 0.10025934231377061, 0.14227750026675529, 0.027846845642827887, 0.06587255617950803, 0.0377115142788618, 0.1709282163040912, 0.1258409006202471, 0.029440758752914818, 0.011341041598241662]
  Cost_jalihal = [0.04352060436354394, 0.16782352005896917, 0.5803849015840222, 0.0813183538847046, 0.02664445093902895, 0.01521022589708537, 0.6459549142857574, 0.6638851878171327, 0.054569323749060414, 0.01492059687602178]

  His_Glucose_addition_Mig1 = His1 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[1], Cost_alt_2[1], Cost_jalihal[1]], 
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "log(Mig1/(Mig1_T-Mig1)",
    legend = false
  )

  His_Glucose_addition_Sch9 = His2 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[2], Cost_alt_2[2], Cost_jalihal[2]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glucose addition Sch9\\Delta",
    legend = false
  )

  His_Glucose_addition_cAMP = His3 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[3], Cost_alt_2[3], Cost_jalihal[3]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glucose addition cAMP",
    legend = false
  )

  His_Glucose_addition_Sch9_p = His4 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[4], Cost_alt_2[4], Cost_jalihal[4]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glucose addition Sch9 p",
    legend = false
  )

  His_Glucose_starvation_Snf1_p = His5 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[5], Cost_alt_2[5], Cost_jalihal[5]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glucose starvation Snf1",
    legend = false
  )

  His_Glucose_starvation_Sch9_p = His6 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[6], Cost_alt_2[6], Cost_jalihal[6]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glucose starvation Sch9",
    legend = false
  )

  His_Sch9P_glutamine_L = His8 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[7], Cost_alt_2[7], Cost_jalihal[7]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glutamine addition, Low Glutamine (0.3)",
    legend = false
  )

  His_Sch9P_glutamine_H = His9 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[8], Cost_alt_2[8], Cost_jalihal[8]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glutamine addition, High glutamine",
    legend = false
  )

  His_Glutamine_addition_Sch9_gtr1Delta = His7 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[9], Cost_alt_2[9], Cost_jalihal[9]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Glutamine addition gtr1\\Delta, Glutamine_{ext} = 1.0",
    legend = false
  )

  His_Rapamycin_treatment = His10 = bar(
    [label_alt_1,label_alt_2,label_jalihal],
    [Cost_alt_1[10], Cost_alt_2[10], Cost_jalihal[10]],  
    bar_width = 1,
    color = [Color_alt_1, Color_alt_2, Color_jalihal],
    ylabel = "Cost",
    title = "Rapamycin treatment, TORC1_T = 0",
    legend = false
  )
  plot1 = plot(His1, His2, His3, His4, His5, His6, His7, His8, His9, His10, layout = (5, 2), legend = false, titlefontsize = 14, guidefontsize = 10, guide_position = :left, margin= 20Plots.mm)
  plot!(plot1,size=(1500,2500))
  #display(plot1)
  return plot1
end

