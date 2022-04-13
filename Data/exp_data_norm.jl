include("exp_data.jl")
include("../Code/Model/ODE_methods.jl")

# För Mig1 gjörs ingen normalisering (logaritmisk skala)
data_Mig1_glucose_relief = raw_data_Mig1_glucose_relief

# data för cAMP använder min och max taget från all data för cAMP
data_sch9Delta_cAMP = Minmaxnorm(raw_data_sch9Delta_cAMP, 0.052, 1.227)
data_cAMP = Minmaxnorm(raw_data_cAMP, 0.052, 1.227)

data_Snf1 = Minmaxnorm(raw_data_Snf1)

data_Sch9_glucose_relief = Minmaxnorm(raw_data_Sch9_glucose_relief) # Rel. fosforylering
data_Sch9_glucose_starve = Minmaxnorm(raw_data_Sch9_glucose_starve) # Rel. fosforylering

# data för Sch9 med enhet %fosforylering har min och max taget från all data för Sch9 för 
data_Sch9_gtr1Delta = Minmaxnorm(raw_data_Sch9_gtr1Delta, 4.82, 52.27)
data_Sch9P_glutamine_L = Minmaxnorm(raw_data_Sch9P_glutamine_L, 4.82, 52.57)
data_Sch9P_glutamine_H = Minmaxnorm(raw_data_Sch9P_glutamine_H, 4.82, 52.57)

data_Rib_rap = Minmaxnorm(raw_data_Rib_rap)
