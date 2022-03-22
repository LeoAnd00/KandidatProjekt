# Rådata för plottning av tidsseriedata. 
# Ev. skalningar av data görs vid plottning (förutom ev. tidsomvandling till min)

# Glucose starvation (Snf1)
t_Snf1=[0.39, 5.55, 30.56, 60.64] # min
data_Snf1=[0.03, 0.26, 0.32, 0.302] # Skalas med maxvärde

# Glucose addition (cAMP)
t_cAMP=[0.16, 15.27, 30.66, 44.88, 59.95, 74.84, 89.72, 104.2, 119.75, 149.39, 180.06]./60 # min
data_cAMP=[0.079, 0.632, 0.822, 1.197, 1.227, 0.928, 0.617, 0.617, 0.632, 0.209, 0.282] # Skalas med maxvärde

# Glucose addition, Sch9delta mutation (cAMP)
t_sch9Delta_cAMP=[0.15, 14.97, 30.09, 44.76, 59.71, 74.74, 89.51, 104.54, 119.59, 149.79, 180.28]./60 # min
data_sch9Delta_cAMP=[0.052, 0.112, 0.67, 0.259, 0.484, 0.375, 0.264, 0.279, 0.301, 0.275, 0.259] # Skalas med maxvärde/2 *känns konstigt - kolla Jalihal's kod

# High glutamine addition (Sch9)
t_Sch9P_glutamine_H=[0.06, 0.36, 0.73, 1.05, 3.07, 4.12, 5.06, 8.05, 11.03, 15.07, 30.06] #min
data_Sch9P_glutamine_H=[5.92, 10.23, 34.29, 41.14, 40.59, 39.82, 32.99, 13.09, 17.25, 35.4, 52.57] # Skalas med maxvärde *Datapunkt vid andra böjen blir lite konstig

# Low glutamine addition (Sch9)
t_Sch9P_glutamine_L=[0.46, 1.09, 1.57, 2.04, 3.01, 4.05, 5.03, 8.05, 11.04, 15.01, 29.98] #min
data_Sch9P_glutamine_L=[10.23, 49.55, 45.84, 27.66, 6.55, 5.49, 5.22, 5.27, 4.82, 4.89, 5.23] # Skalas med maxvärde *Blir inte rätt. Skalas med maximal modellvärde? Hur?

# Glucose addition, gtr1Delta mutation (Sch9)
t_Sch9_gtr1Delta=[0.15, 0.7, 1.04, 1.51, 2.06, 3.05, 4.04, 5.06, 8.06, 11.0, 15.04, 30.02] #min
data_Sch9_gtr1Delta=[7.96, 8.48, 8.48, 7.44, 7.44, 7.31, 7.44, 7.96, 9.64, 12.49, 25.7, 51.72] # Skalas med maxvärde *Blir inte rätt. Första värdena behöver neråt

# Glucose starvation (Sch9)
t_Sch9_glucose_starve=[0.0,2.5,5.0,15.0,30.0] #min
data_Sch9_glucose_starve=[1.0,0.34,0.13,0.15,0.20] # Hur skalar man denna?

# Glucose relief (Sch9)
t_Sch9_glucose_relief=[0.0,2.0,5.0,15.0,30.0] # min
data_Sch9_glucose_relief=[0.20,0.63,0.92,0.91,1.0] # Skalas ej. *Första punkten är jättekonstig

# Glucose relief (log(nMig1/cMig1))
t_Mig1_glucose_relief=[0.0, 67.0, 104.0, 130.0, 233.0, 334.0, 466.0, 606.0, 771.0, 919.0, 1088.0]./60 # min
data__Mig1_glucose_relief=[1.070, 1.233, 1.363, 1.428, 1.490, 1.465, 1.432, 1.425, 1.416, 1.400, 1.402] # Ingen skalning. 

# Rapamycin treament (Rel. Rib)
t_Rib_rap=[0.0,15.0,30.0,60.0,90.0] #min
data_Rib_rap=[1.00,0.43,0.19,0.09,0.07] # Ingen skalning. De flesta punkterna ser bra ut

# t_Gln3_rap=[0,5,10,15,20] #min
# data_Gln3_rap=[0.0, 1048.464, 7846.120999999999, 9644.262999999999, 12018.475999999999]