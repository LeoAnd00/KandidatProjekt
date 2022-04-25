# Creates the simulations from the Perturbation-Experiment for both wild typ and mutant
include("perturbation_methods.jl")

# Perturbation simulations Snf1
Snf1_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_glucose_starvation]), last.(nutrient_shifts[index_glucose_starvation]), Snf1)
Snf1_mutant_shift = Steady_state_pertubation_solver([Sak_T => 0.0], first.(nutrient_shifts[index_glucose_starvation]), last.(nutrient_shifts[index_glucose_starvation]), Snf1)

# Gis1
Gis1_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_glucose_addition]), last.(nutrient_shifts[index_glucose_addition]), Gis1)
Gis1_mutant_shift = Steady_state_pertubation_solver([Sch9_T => 0.0], first.(nutrient_shifts[index_glucose_addition]), last.(nutrient_shifts[index_glucose_addition]), Gis1)
# Nth1 (Trehalas)
Tre_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_glucose_starvation]), last.(nutrient_shifts[index_glucose_starvation]), Trehalase)
Tre_mutant_shift = Steady_state_pertubation_solver([PKA_T => 0.3], first.(nutrient_shifts[index_glucose_starvation]), last.(nutrient_shifts[index_glucose_starvation]), Trehalase)

# cAMP
cAMP_shift_pde = Steady_state_pertubation_solver(first.(nutrient_shifts[index_glucose_addition]), last.(nutrient_shifts[index_glucose_addition]), cAMP)
cAMP_mutant_pde_shift = Steady_state_pertubation_solver([PDE_T => 0.0], first.(nutrient_shifts[index_glucose_addition]), last.(nutrient_shifts[index_glucose_addition]), cAMP)	
cAMP_shift_ras = Steady_state_pertubation_solver(first.(nutrient_shifts[index_glucose_addition]), last.(nutrient_shifts[index_glucose_addition]), cAMP)
cAMP_mutant_ras_shift = Steady_state_pertubation_solver([Ras_T => 0.0, w_pka_camp => 0.0], first.(nutrient_shifts[index_glucose_addition]), last.(nutrient_shifts[index_glucose_addition]), cAMP)

# Rib
Rib_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Rib)
Rib_mutant_shift = Steady_state_pertubation_solver([Sch9_T => 0.0], first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Rib)

# Sch9 wild typ
Sch9_lst4_lst7_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Sch9)
Sch9_gtr1_gtr2_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Sch9)
Sch9_gtr1_2_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Sch9)

# Sch9 mutant
Sch9_lst4_lst7_mutant_shift = Steady_state_pertubation_solver([EGOGAP_T => 0.0], first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Sch9)
Sch9_gtr1_gtr2_mutant_shift = Steady_state_pertubation_solver([EGO_T => 0.0, w_torc_ego => 0.0, w_torc_egoin => 0.0, w_torc_glut => 0.5], first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Sch9)
Sch9_gtr1_2_mutant_shift = Steady_state_pertubation_solver([EGO_T => 0.0, w_torc_ego => 0.0, w_torc_egoin => 0.0], first.(nutrient_shifts[index_high_glutamine]), last.(nutrient_shifts[index_high_glutamine]), Sch9)
 
# Gcn4
Gcn4_shift = Steady_state_pertubation_solver(first.(nutrient_shifts[index_nitrogen_starvation]), last.(nutrient_shifts[index_nitrogen_starvation]), Gcn4)
Gcn4_mutant_shift = Steady_state_pertubation_solver([Gcn2_T => 0.0], first.(nutrient_shifts[index_nitrogen_starvation]), last.(nutrient_shifts[index_nitrogen_starvation]), Gcn4)

# Gln3
Gln3_shift_sit = Steady_state_pertubation_solver(first.(nutrient_shifts[index_nitrogen_starvation]), last.(nutrient_shifts[index_nitrogen_starvation]), Gln3)
Gln3_mutant_shift_sit = Steady_state_pertubation_solver([w_gln_sit => 0.0, TORC1_T => 0.0], first.(nutrient_shifts[index_nitrogen_starvation]), last.(nutrient_shifts[index_nitrogen_starvation]), Gln3)

# Normalisation of Rib and cAMP simulations
cAMP_shifts = [cAMP_shift_ras, cAMP_shift_pde, cAMP_mutant_ras_shift, cAMP_mutant_pde_shift]
min_cAMP = minimum([first.(cAMP_shifts) last.(cAMP_shifts)])
max_cAMP = maximum([first.(cAMP_shifts) last.(cAMP_shifts)])

Rib_shifts = [Rib_shift, Rib_mutant_shift]
min_Rib = minimum([first.(Rib_shifts) last.(Rib_shifts)])
max_Rib = maximum([first.(Rib_shifts) last.(Rib_shifts)])


cAMP_shifts_wt = [cAMP_shift_pde, cAMP_shift_ras]
cAMP_shifts_wt_norm = Minmaxnorm_shifts(cAMP_shifts_wt, min_cAMP, max_cAMP)

Rib_shifts_wt = [Rib_shift]
Rib_shifts_wt_norm = Minmaxnorm_shifts(Rib_shifts_wt, min_Rib, max_Rib)

raw_shifts_wt = [Snf1_shift, Gis1_shift, Tre_shift, Sch9_gtr1_2_shift, Gcn4_shift, Gln3_shift_sit, 
                 Sch9_gtr1_gtr2_shift, Sch9_lst4_lst7_shift
]

# Normalisation of the expeimental data
# The first pair is wild typ and the next pair is mutant for each simulation
exp_data_Rib = [150.0 => 600.0, 100.0 => 100.0]
exp_data_Rib_norm = Minmaxnorm_data(exp_data_Rib)

exp_data_cAMP_pde = [0.25 => 1.1, 0.25 => 3.8]
exp_data_cAMP_pde_norm = Minmaxnorm_data(exp_data_cAMP_pde)

exp_data_cAMP_ras = [0.25 => 1.1, 0.25 => 0.25]
exp_data_cAMP_ras_norm = Minmaxnorm_data(exp_data_cAMP_ras, 0.0, 3.8)

exp_data_Snf1 = [0.05 => 2.25, 0.01 => 0.01]
exp_data_Snf1_norm = Minmaxnorm_data(exp_data_Snf1)

exp_data_Gis1 = [700.0 => 35.0, 100.0 => 10.0]
exp_data_Gis1_norm = Minmaxnorm_data(exp_data_Gis1)

exp_data_Tre = [100.0 => 25.0, 50.0 => 25.0]
exp_data_Tre_norm = Minmaxnorm_data(exp_data_Tre)

exp_data_Sch9_gtr1_2 = [5.0 => 55.0, 5.0 => 55.0]
exp_data_Sch9_gtr1_2_norm = Minmaxnorm_data(exp_data_Sch9_gtr1_2)

exp_data_Gcn4 = [7.3 => 100.0, 11.0 => 38.0]
exp_data_Gcn4_norm = Minmaxnorm_data(exp_data_Gcn4)

exp_data_Gln3_sit = [461.0 => 6248.0, 0.0 => 789.0]
exp_data_Gln3_sit_norm = Minmaxnorm_data(exp_data_Gln3_sit)

exp_data_wt_norm = vcat(exp_data_Rib_norm[1], exp_data_cAMP_pde_norm[1], exp_data_cAMP_ras_norm[1], exp_data_Snf1_norm[1], exp_data_Gis1_norm[1], exp_data_Tre_norm[1], exp_data_Sch9_gtr1_2_norm[1], exp_data_Gcn4_norm[1], exp_data_Gln3_sit_norm[1])

# Plotta resultatet in a bar-plott
using Plots
using StatsPlots

######## Wild typ ################
shifts_wt = vcat(Rib_shifts_wt_norm, cAMP_shifts_wt_norm, raw_shifts_wt)

# shift_names_wt = ["\$\\mathrm{Rib}\$", "\$\\mathrm{cAMP_{pde}}\$", "\$\\mathrm{cAMP_{ras}}\$", "\$\\mathrm{Snf1}\$", "\$\\mathrm{Gis1}\$", "\$\\mathrm{Nth1}\$", "\$\\mathrm{Sch9_{gtr1-2}}\$","\$\\mathrm{Gcn4}\$", "\$\ \mathrm{Gln3_{sit}}\$", "\$\\mathrm{Sch9_{gtr1-gtr2}}\$", "\$\\mathrm{Sch9_{lst4-1st7}}\$"]
shift_names_wt = ["\$\\mathrm{Rib}\$", "\$\\mathrm{cAMP}\$", "\$\\mathrm{cAMP}\$", 
 "\$\\mathrm{Snf1}\$", "\$\\mathrm{Gis1}\$", "\$\\mathrm{Nth1}\$", 
 "\$\\mathrm{Sch9_{gtr1-2}}\$", "\$\\mathrm{Gcn4}\$", "\$\\mathrm{Gln3}\$", 
 "\$\\mathrm{Sch9_{gtr1-gtr2}}\$","\$\\mathrm{Sch9_{lst4-lst7}}\$"
]

# Calculates the diffrence between steady state solutions, add them to an epty vector
shift_magnitude_wt = AbstractFloat[]
for i in range(1, length(shifts_wt))
    push!(shift_magnitude_wt, last.(shifts_wt[i]) - first.(shifts_wt[i]))
end 
 
# Instead of using the highest and lowest value from the data and simulation, 0 and 100 is used as scaling factors 
exp_data_wt_special_raw = [2.0 => 73.0, 2.0 => 73.0] # Sch9_gtr1_gtr2, Sch9_lst4_lst7
exp_data_wt_special_norm = Minmaxnorm_shifts(exp_data_wt_special_raw, 0, 100.0)

exp_data_wt_raw = [150.0 => 600.0, 0.25 => 1.1, 0.25 => 1.1, 0.05 => 2.25, 700.0 => 35.0, 100.0 => 25.0, 5.0 => 55.0, 7.3 => 100.0, 461.0 => 6248.0]  

# Combinde all experimental data
exp_data_wt = vcat(exp_data_wt_norm, exp_data_wt_special_norm)

 # Calculates the diffrence between experimental data
exp_shift_magnitude_wt = AbstractFloat[]
for i in range(1, length(exp_data_wt))
    push!(exp_shift_magnitude_wt, last.(exp_data_wt[i]) - first.(exp_data_wt[i]))
end

# Creats a vector with the diffrence from pre- and post-shift for both experimental data and simulations ready to plot
shift_values_wt = [shift_magnitude_wt exp_shift_magnitude_wt] 

plot1 = groupedbar(                                                                                                                                         
    shift_values_wt, group=repeat(["Simulering", "Experiment"], inner = length(shift_names_wt)), xlims=(-1.05, 1.05),
    yticks=(1:length(shift_names_wt), shift_names_wt), bar_position=:group, ytickfont=font(9), legend=:topleft,
    ylabel="Readouts, wt", xlabel="Shift magnitude", bar_width = 0.65, framestyle=:box, orientation=:horizontal
)
display(plot1)

####### MUTANT #######
# DO the same thing with de mutants

# Scales the data according to the Minmaxnorm function 
# cAMP
cAMP_shifts_mutant = [cAMP_mutant_pde_shift, cAMP_mutant_ras_shift]
cAMP_shifts_mutant_norm = Minmaxnorm_shifts(cAMP_shifts_mutant, min_cAMP, max_cAMP)

# Rib
Rib_shifts_mutant = [Rib_mutant_shift]
Rib_shifts_mutant_norm = Minmaxnorm_shifts(Rib_shifts_mutant, min_Rib, max_Rib)

# Placing the mutant simulations i order for the plot
raw_shifts_mutant = [Snf1_mutant_shift, Gis1_mutant_shift, Tre_mutant_shift, Sch9_gtr1_2_mutant_shift, Gcn4_mutant_shift, Gln3_mutant_shift_sit,
Sch9_gtr1_gtr2_mutant_shift, Sch9_lst4_lst7_mutant_shift]
 
shifts_mutant = vcat(Rib_shifts_mutant_norm, cAMP_shifts_mutant_norm, raw_shifts_mutant)

shift_names_mutant = ["\$\\mathrm{Rib_{sch9}}\$", "\$\\mathrm{cAMP_{pde1}}\$", "\$\\mathrm{cAMP_{ras1ras2bcy1}}\$", 
 "\$\\mathrm{Snf1_{sak1tos3elm1}}\$", "\$\\mathrm{Gis1_{sch9}}\$", "\$\\mathrm{Nth1_{tpk1tpk2}}\$", 
 "\$\\mathrm{Sch9_{gtr1-2}}\$", "\$\\mathrm{Gcn4_{gcn2}}\$", "\$\\mathrm{Gln3_{sit}}\$", 
 "\$\\mathrm{Sch9_{gtr1-gtr2}}\$","\$\\mathrm{Sch9_{lst4-lst7}}\$"]

# Calculates the diffrence between steady state solutions, add them to an epty vector
shift_magnitude_mutant = AbstractFloat[]
for i in range(1, length(shifts_mutant))
    push!(shift_magnitude_mutant, last.(shifts_mutant[i]) - first.(shifts_mutant[i]))
end

# Instead of using the highest and lowest value from the data and simulation, 0 and 100 is used as scaling factors 
exp_data_mutant_special_raw = [22.0 => 12.0, 1.0 => 19.0] # Sch9_gtr1_gtr2_mutant, # Sch9_lst4_lst7_mutant
exp_data_mutant_special_norm = Minmaxnorm_shifts(exp_data_mutant_special_raw, 0.0, 100.0)

exp_data_mutant_raw = [100.0 => 100.0, 0.25 => 3.8, 0.25 => 0.25, 0.01 => 0.01, 100.0 => 10.0, 
50.0 => 25.0, 5.0 => 55.0, 11.0 => 38.0, 0.0 => 789.0] 

exp_data_mutant_norm = vcat(exp_data_Rib_norm[2], exp_data_cAMP_pde_norm[2], exp_data_cAMP_ras_norm[2], exp_data_Snf1_norm[2], exp_data_Gis1_norm[2], exp_data_Tre_norm[2], exp_data_Sch9_gtr1_2_norm[2], exp_data_Gcn4_norm[2], exp_data_Gln3_sit_norm[2])

exp_data_mutant = vcat(exp_data_mutant_norm, exp_data_mutant_special_norm)

# Calculates the diffrence between steady state solutions, add them to an epty vector
exp_shift_magnitude_mutant = AbstractFloat[]
for i in range(1, length(exp_data_mutant))
    push!(exp_shift_magnitude_mutant, last.(exp_data_mutant[i]) - first.(exp_data_mutant[i]))
end

# Creats the final vector, then plot results
shift_values_mutant = [shift_magnitude_mutant exp_shift_magnitude_mutant]

plot2 = groupedbar(
    shift_values_mutant, group=repeat(["Simulering", "Experiment"], inner = length(shift_names_mutant)),
    yticks=(1:length(shift_names_mutant), shift_names_mutant), xlims=(-1.05, 1.05), bar_position=:group, bar_width=0.65, legend=:topleft,
    ylabel="Readouts mutanter", xlabel="Shift magnitude", orientation=:horizontal, ytickfont = font(9), framestyle=:box
)
display(plot2)