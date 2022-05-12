include("simulation_methods.jl")

default(titlefontsize=5)
default(xtickfont=4)
default(ytickfont=4)
default(guidefont=6)
default(ylabelfontsize=5)
default(xlabelfontsize=5)
default(lw=0.6)

###### Glucose starvation ######
sol_array_gluc_starve = Generate_sol_array(index_glucose_starvation, (0.0, 40.0)) # Generates array of simulation

p1=return_simulation_result(sol_array_gluc_starve, Mig1, (0.0, 20.0), (0.8, 1.0), "Glukossvältning", :topright) 
p2=return_simulation_result(sol_array_gluc_starve, Snf1, (0.0, 40.0), (0.0, 1.1), "Glukossvältning", :right) 
p3=return_simulation_result(sol_array_gluc_starve, Gis1, (0.0, 10.0), (0.0, 1.0), "Glukossvältning", :topright) 
p4=return_simulation_result(sol_array_gluc_starve, Gln3, (0.0, 40.0), (0.0, 1.1), "Glukossvältning", :right) 
p5=return_simulation_result(sol_array_gluc_starve, Dot6, (0.0, 4.0), (0.9, 1.0), "Glukossvältning",  :topright) 


###### Nitrogen starvation ######
sol_array_nit_starve = Generate_sol_array(index_nitrogen_starvation, (0.0, 40.0))

p6=return_simulation_result(sol_array_nit_starve, Snf1, (0.0, 15.0), (0.0, 0.10), "Kvävesvältning", :right) 
p7=return_simulation_result(sol_array_nit_starve, TORC1, (0.0, 40.0), (0.0, 1.1), "Kvävesvältning", :right) 
p8=return_simulation_result(sol_array_nit_starve, Sch9, (0.0, 40.0), (0.0, 1.1), "Kvävesvältning", :right) 
p9=return_simulation_result(sol_array_nit_starve, Gln3, (0.0, 40.0), (0.0, 1.1), "Kvävesvältning", :right) 
p10=return_simulation_result(sol_array_nit_starve, Rtg13, (0.0, 40.0), (0.0, 1.1), "Kvävesvältning", :bottomright)
p11=return_simulation_result(sol_array_nit_starve, Dot6, (0.0, 10.0), (0.9, 1.0), "Kvävesvältning",:bottomright) 


###### Glucose addition ######
sol_array_gluc_add = Generate_sol_array(index_glucose_addition, (0.0, 15.0))

p12=return_simulation_result(sol_array_gluc_add, Snf1, (0.0, 10.0), (0.0, 1.0), "Glukostillsättning", :right) 
p13=return_simulation_result(sol_array_gluc_add, Gis1, (0.0, 7.0), (0.0, 1.0), "Glukostillsättning", :right) 

###### Gtr1Delta | EGO deletion, Nitrogen addition ######
sol_array_nit_add = Generate_sol_array(index_high_glutamine, (0.0, 40.0))
sol_array_gtr1Delta = Generate_sol_array(index_high_glutamine, (0.0, 40.0), [EGO_T => 0.0])

p14=return_wt_mut_results(
    sol_array_nit_add, sol_array_gtr1Delta, Gln3, (0.0, 25.0), (0.0, 1.0), 
    "Kvävetillsättning, gtr1\$\\Delta\$", :right
    )

###### Snf1 deletion, Glucose starvation ######
sol_array_Snf1Delta = Generate_sol_array(index_glucose_starvation, (0.0, 1.0), [Snf1_T => 0.0])

p15=return_wt_mut_results(
    sol_array_gluc_starve, sol_array_Snf1Delta, Gln3, (0.0, 0.5), (0.0, 0.15), 
    "Glukossvältning, Snf1\$\\Delta\$", :right
    )
p16=return_wt_mut_results(
    sol_array_gluc_starve, sol_array_Snf1Delta, Cyr1, (0.0, 0.5), (0.0, 0.7), 
    "Glukossvältning, Snf1\$\\Delta\$", :right
    )

plot_list = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16]
l = @layout [
    grid(4,4)
]

plot_all = plot(plot_list..., layout=l, legend = false)
display(plot_all)

file_name = "model_analysis_subplot"
savefig(plot_all, pwd()*"/Results/Model_analysis/"*file_name*".pdf")
savefig(plot_all, pwd()*"/Results/Model_analysis/"*file_name*".png")
