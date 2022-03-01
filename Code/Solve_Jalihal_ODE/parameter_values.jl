# Definierar optimala parametervärden (funna i Table S2).
# Ny rad motsvarar parametrar för ny ODE enligt kommentarer
# För model inputs är ATP, Carbon (Glukos) och Glutamine_ext = 1. NH4 och Proline är = 0
# "Optimala" inputs är dock då Carbon, ATP och Glutamine_ext är 1
# gamma_gcn4, gamma_rtg13, gamma_gis1, gamma_dot6 antas vara ett.

p = [gammacyr => 8.96, Cyr1_T => 1.0, w_cyr_glu => 5.13, w_cyr_snf => 0.12, sigma_cyr => 3.5, w_cyr => 1.35, # Cyr1
w_dot => 0.29, sigma_dot => 20.0, w_dot_sch_pka => 0.16, Dot6_T => 1.0, # Dot6
w_ego => 0.28, sigma_ego => 5.0, w_ego_basal => 0.01, EGO_T => 1.0, gammaego => 50.66, w_ego_gap => 2.21, # EGO
gammagap => 0.56, sigma_gap => 1.0, w_gap_torc => 88.33, w_gap_N => 7.76, EGOGAP_T => 1.0, # EGOGAP
gamma_gcn2 => 4.71, Gcn2_T => 1.0, w_gcn_torc => 1.29, sigma_gcn2 => 20.0, w_gcn => 0.12, # Gcn2
w_gcn4_gcn2_trna => 1.53, w_gcn4 => 0.74, tRNA_sensitivity => 74.51, sigma_gcn4 => 5.0, tRNA_total => 2.47, Gcn4_T => 1.0, # Gcn4
sigma_gis1 => 10.0, w_gis => 1.3, Gis1_T => 1.0, w_gis_pka => 3.3, w_gis_sch => 0.84, #Gis1
w_gln1_gln3 => 0.52, Gln1_T => 1.0, w_gln1 => 0.22, gammagln1 => 0.06, sigma_gln1 => 1.0, # Gln1
w_gln_snf => 3.9, sigma_gln => 10.0, w_gln_sit => 0.86, Gln3_T => 1.0, w_gln3 => 0.64, gammagln3 => 0.08, # Gln3
k_acc_pro => 0.0, k_acc_glu => 0.05, k_degr => 0.09, k_acc_nh4 => 0.0, Glutamine_ext => 1.0, NH4 => 0.0, Carbon => 1.0, Proline => 0.0, # Glutamine
w_mig_snf => 1.21, w_mig => 10.64, sigma_mig1 => 0.27, Mig1_T => 1.0, gamma_mig => 0.66, w_mig_pka => 2.31, # Mig1
sigma_pde => 1.9, gammapde => 0.28, w_pde_pka => 2.89, PDE_T => 1.0, w_pde => 0.38, # PDE
w_pka => 0.06, w_pka_sch9 => 17.5, PKA_T => 1, w_pka_camp => 102.11, sigma_pka => 1.0, gammapka => 2.68, # PKA
k_pr => 0.02, ATP => 1.0, # Protein
w_ras_pka => 1.87, sigma_ras => 1.0, Ras_T => 1.0, gammaras => 1.82, w_ras_glu => 0.21, w_ras => 0.02, # Ras
k_mRNA_degr => 0.07, k_transcription => 0.24, # Rib
w_rtg => 0.19, sigma_rtg => 10.0, w_rtg_torc => 0.88, Rtg13_T => 1.0, # Rtg13
w_sak => 0.21, Sak_T => 1.0, w_sak_pka => 0.38, sigma_sak => 20.0, # Sak
w_sch9_torc => 1.96, w_sch9 => 0.57, gammasch9 => 4.63, sigma_sch9 => 8.0, Sch9_T => 1.0, # Sch9
w_snf_sak => 1.52, gammasnf => 0.82, w_snf_glc => 1.15, sigma_snf => 3.0, w_snf => 0.54, Snf1_T => 1.0, # Snf1
w_torc_egoin => 0.3, TORC1_T => 1.0, w_torc_ego => 0.88, w_torc => 0.54, w_torc_snf => 0.44, w_torc_glut => 0.86, sigma_tor => 5.0, gammator => 7.55, # TORC1
w_tps_pka => 0.57, sigma_tps => 5.0, gammatps => 0.47, w_tps => 0.05, PKA_T => 1.0, Tps1_T => 1.0, # Tps1
w_tre => 1.07, gammatre => 0.34, w_tre_pka => 3.07, Trehalase_T => 1.0, sigma_trehalase => 10.0, # Trehalase
k_camp_cyr => 10.87, k_camp_deg => 0.08, k_camp_pde => 14.12, # cAMP
eIF_T => 1.0, w_eif_gcn2 => 0.28, sigma_eif => 1.0, gammaeif => 0.47, w_eif => 3.73] # eIF

# Look up table, gamma_gcn4, gamma_rtg13, gamma_gis1, gamma_dot6 saknas i denna modellen

p_lookup_table = ["gammacyr", "Cyr1_T", "w_cyr_glu", "w_cyr_snf", "sigma_cyr", "w_cyr", # Cyr1
"w_dot", "sigma_dot", "w_dot_sch_pka", "Dot6_T", # Dot6
"w_ego", "sigma_ego", "w_ego_basal", "EGO_T", "gammaego", "w_ego_gap", # EGO
"gammagap", "sigma_gap", "w_gap_torc", "w_gap_N", "EGOGAP_T", # EGOGAP
"gamma_gcn2", "Gcn2_T", "w_gcn_torc", "sigma_gcn2", "w_gcn", # Gcn2
"w_gcn4_gcn2_trna", "w_gcn4", "tRNA_sensitivity", "sigma_gcn4", "tRNA_total", "Gcn4_T", # Gcn4
"sigma_gis1", "w_gis", "Gis1_T", "w_gis_pka", "w_gis_sch", #Gis1
"w_gln1_gln3", "Gln1_T", "w_gln1", "gammagln1", "sigma_gln1", # Gln1
"w_gln_snf", "sigma_gln", "w_gln_sit", "Gln3_T", "w_gln3", "gammagln3", # Gln3
"k_acc_pro", "k_acc_glu", "k_degr", "k_acc_nh4", "Glutamine_ext", "NH4", "Carbon", "Proline", # Glutamine
"w_mig_snf", "w_mig", "sigma_mig1", "Mig1_T", "gamma_mig", "w_mig_pka", # Mig1
"sigma_pde", "gammapde", "w_pde_pka", "PDE_T", "w_pde", # PDE
"w_pka", "w_pka_sch9", "PKA_T", "w_pka_camp", "sigma_pka", "gammapka", # PKA
"k_pr", "ATP", # Protein
"w_ras_pka", "sigma_ras", "Ras_T", "gammaras", "w_ras_glu", "w_ras", # Ras
"k_mRNA_degr", "k_transcription", # Rib
"w_rtg", "sigma_rtg", "w_rtg_torc", "Rtg13_T", # Rtg13
"w_sak", "Sak_T", "w_sak_pka", "sigma_sak", # Sak
"w_sch9_torc", "w_sch9", "gammasch9", "sigma_sch9", "Sch9_T", # Sch19
"w_snf_sak", "gammasnf", "w_snf_glc", "sigma_snf", "w_snf", "Snf1_T", # Snf1
"w_torc_egoin", "TORC1_T", "w_torc_ego", "w_torc", "w_torc_snf", "w_torc_glut", "sigma_tor", "gammator", # TORC1
"w_tps_pka", "sigma_tps", "gammatps", "w_tps", "PKA_T", "Tps1_T", # Tps1
"w_tre", "gammatre", "w_tre_pka", "Trehalase_T", "sigma_trehalase", # Trehalase
"k_camp_cyr", "ATP", "k_camp_deg", "k_camp_pde", # cAMP
"eIF_T", "w_eif_gcn2", "sigma_eif", "gammaeif", "w_eif"] # eIF