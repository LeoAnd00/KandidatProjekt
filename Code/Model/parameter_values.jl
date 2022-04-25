 # Defines optimal parametervalues from Jalihal et. all
# A new row reprecent parameters for the new ODE according to the comments below
#  The model imputs ATP, Carbon, Glutamine_ext are all equal to 1, NH4 and PRoline are equal to 0
# The optimal imputs are when Carbon, ATP and Glutamine_ext are equal to 1


p_var = [gammacyr => 8.96, w_cyr_glu => 5.13, w_cyr_snf => 0.12, w_cyr => 1.35, # Cyr1
w_dot => 0.29, w_dot_sch_pka => 0.16, # Dot6
w_ego => 0.28, w_ego_basal => 0.01, gammaego => 50.66, w_ego_gap => 2.21, # EGO
gammagap => 0.56, w_gap_torc => 88.33, w_gap_N => 7.76, # EGOGAP
gamma_gcn2 => 4.71, w_gcn_torc => 1.29, w_gcn => 0.12, # Gcn2
w_gcn4_gcn2_trna => 1.53, w_gcn4 => 0.74, tRNA_sensitivity => 74.51, tRNA_total => 2.47, # Gcn4
w_gis => 1.3, w_gis_pka => 3.3, w_gis_sch => 0.84, #Gis1
w_gln1_gln3 => 0.52, w_gln1 => 0.22, gammagln1 => 0.06, # Gln1
w_gln_snf => 3.9, w_gln_sit => 0.86, w_gln3 => 0.64, gammagln3 => 0.08, # Gln3
k_acc_pro => 0.0, k_acc_glu => 0.05, k_degr => 0.09, k_acc_nh4 => 0.0, # Glutamine
w_mig_snf => 1.21, w_mig => 10.64, gamma_mig => 0.66, w_mig_pka => 2.31, # Mig1
gammapde => 0.28, w_pde_pka => 2.89, w_pde => 0.38, # PDE
w_pka => 0.06, w_pka_sch9 => 17.5, w_pka_camp => 102.11, gammapka => 2.68, # PKA
k_pr => 0.02, # Protein
w_ras_pka => 1.87, gammaras => 1.82, w_ras_glu => 0.21, w_ras => 0.02, # Ras
k_mRNA_degr => 0.07, k_transcription => 0.24, # Rib
w_rtg => 0.19, w_rtg_torc => 0.88, # Rtg13
w_sak => 0.21, w_sak_pka => 0.38, # Sak
w_sch9_torc => 1.96, w_sch9 => 0.57, gammasch9 => 4.63, # Sch9
w_snf_sak => 1.52, gammasnf => 0.82, w_snf_glc => 1.15, w_snf => 0.54, # Snf1
w_torc_egoin => 0.3, w_torc_ego => 0.88, w_torc => 0.54, w_torc_snf => 0.44, w_torc_glut => 0.86, gammator => 7.55, # TORC1
w_tps_pka => 0.57, gammatps => 0.47, w_tps => 0.05, # Tps1
w_tre => 1.07, gammatre => 0.34, w_tre_pka => 3.07, # Trehalase
k_camp_cyr => 10.87, k_camp_deg => 0.08, k_camp_pde => 14.12, # cAMP
w_eif_gcn2 => 0.28, gammaeif => 0.47, w_eif => 3.73] # eIF

p_var_names = first.(p_var)

# Create look-up table
p_var_lookup_table = []
for i in range(1, length(p_var_names))
    push!(p_var_lookup_table, string(p_var_names[i]))
end

p_const = [Cyr1_T => 1.0, sigma_cyr => 3.5, # Cyr1
sigma_dot => 20.0, Dot6_T => 1.0, # Dot6
sigma_ego => 5.0, EGO_T => 1.0, # EGO
sigma_gap => 1.0, EGOGAP_T => 1.0, # EGOGAP
Gcn2_T => 1.0, sigma_gcn2 => 20.0, # Gcn2
sigma_gcn4 => 5.0, Gcn4_T => 1.0, # Gcn4
sigma_gis1 => 10.0, Gis1_T => 1.0, #Gis1
Gln1_T => 1.0, sigma_gln1 => 1.0, # Gln1
sigma_gln => 10.0, Gln3_T => 1.0, # Gln3
NH4 => 0.0, Proline => 0.0, # Glutamine
sigma_mig1 => 0.27, Mig1_T => 1.0, # Mig1
sigma_pde => 1.9, PDE_T => 1.0, # PDE
PKA_T => 1.0, sigma_pka => 1.0, # PKA
# Protein
sigma_ras => 1.0, Ras_T => 1.0, # Ras
# Rib
sigma_rtg => 10.0, Rtg13_T => 1.0, # Rtg13
Sak_T => 1.0, sigma_sak => 20.0, # Sak
sigma_sch9 => 8.0, Sch9_T => 1.0, # Sch9
sigma_snf => 3.0, Snf1_T => 1.0, # Snf1
TORC1_T => 1.0, sigma_tor => 5.0, # TORC1
sigma_tps => 5.0, Tps1_T => 1.0, # Tps1
Trehalase_T => 1.0, sigma_trehalase => 10.0, # Trehalase
# cAMP
eIF_T => 1.0, sigma_eif => 1.0] # eIF

p_const_names = first.(p_const)

p_const_lookup_table = []
for i in range(1, length(p_const_names))
    push!(p_const_lookup_table, string(p_const_names[i]))
end

p_conc = vcat(p_const, p_var)
p_conc_lookup_table = vcat(p_const_lookup_table, p_var_lookup_table)