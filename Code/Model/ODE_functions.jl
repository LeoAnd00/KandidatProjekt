# Defines the ODE-model in ModelingToolkit

using ModelingToolkit
using DifferentialEquations

# Defines parameters, variables and derivative operator for constructing the model using ModelingToolkit
@parameters gammacyr Cyr1_T w_cyr_glu w_cyr_snf sigma_cyr w_cyr w_dot sigma_dot w_dot_sch_pka Dot6_T w_ego sigma_ego w_ego_basal EGO_T gammaego w_ego_gap gammagap sigma_gap w_gap_torc w_gap_N EGOGAP_T gamma_gcn2 Gcn2_T w_gcn_torc sigma_gcn2 w_gcn w_gcn4_gcn2_trna w_gcn4 tRNA_sensitivity sigma_gcn4 tRNA_total Gcn4_T sigma_gis1 w_gis Gis1_T w_gis_pka w_gis_sch w_gln1_gln3 Gln1_T w_gln1 gammagln1 sigma_gln1 w_gln_snf sigma_gln w_gln_sit Gln3_T w_gln3 gammagln3 k_acc_pro k_acc_glu k_degr k_acc_nh4 NH4 Proline w_mig_snf w_mig sigma_mig1 Mig1_T gamma_mig w_mig_pka sigma_pde gammapde w_pde_pka PDE_T w_pde w_pka w_pka_sch9 w_pka_camp sigma_pka gammapka k_pr tRNA_total w_ras_pka sigma_ras Ras_T gammaras w_ras_glu w_ras k_mRNA_degr k_transcription w_rtg sigma_rtg w_rtg_torc Rtg13_T w_sak Sak_T w_sak_pka sigma_sak w_sch9_torc w_sch9 gammasch9 sigma_sch9 Sch9_T w_snf_sak gammasnf w_snf_glc sigma_snf w_snf Snf1_T w_torc_egoin w_torc w_torc_snf w_torc_glut sigma_tor gammator TORC1_T w_torc_ego w_tps_pka sigma_tps gammatps w_tps PKA_T Tps1_T w_tre gammatre w_tre_pka Trehalase_T sigma_trehalase k_camp_cyr ATP k_camp_deg k_camp_pde eIF_T w_eif_gcn2 sigma_eif gammaeif w_eif
@variables t Cyr1(t) Dot6(t) EGO(t) EGOGAP(t) Gcn2(t) Gcn4(t) Gis1(t) Gln1(t) Gln3(t) Glutamine(t) Mig1(t) PDE(t) PKA(t) Protein(t) Ras(t) Rib(t) Rtg13(t) Sak(t) Sch9(t) Snf1(t) TORC1(t) Tps1(t) Trehalase(t) cAMP(t) eIF(t) Carbon(t) ATP(t) Glutamine_ext(t)
D = Differential(t)

function soft_Heaviside(sigma, W)
    1/(1+exp(-sigma*W))
end

# Defines the ODE-system for the implemented model
eqs = [
    # Nutrient signal sensing and transduction 
    D(Glutamine) ~ (k_acc_glu*Glutamine_ext+k_acc_pro*Proline+k_acc_nh4*NH4*Gln1*Carbon)-k_degr*Glutamine,
    D(Cyr1) ~ gammacyr*(Cyr1_T*soft_Heaviside(sigma_cyr, w_cyr_glu*Carbon*Ras-w_cyr-w_cyr_snf*Snf1)-Cyr1), 
    D(Ras) ~ gammaras* (Ras_T*soft_Heaviside(sigma_ras,-w_ras_pka*PKA+w_ras_glu*Carbon+w_ras)-Ras),
    D(EGO) ~ gammaego*(EGO_T*soft_Heaviside(sigma_ego,w_ego_gap*EGOGAP*(Glutamine_ext+0.5*NH4+0.01*Proline)-w_ego*(1-Glutamine)-w_ego_basal)-EGO),
    D(EGOGAP) ~ gammagap*(EGOGAP_T*soft_Heaviside(sigma_gap, w_gap_N*(1-Glutamine)-w_gap_torc*TORC1)-EGOGAP),
    D(cAMP) ~ (k_camp_cyr*Cyr1*ATP-k_camp_pde*PDE*cAMP-k_camp_deg*cAMP),
    D(PDE) ~ gammapde*(PDE_T*soft_Heaviside(sigma_pde,w_pde_pka*PKA-w_pde)-PDE),
    D(Sak) ~ (Sak_T*soft_Heaviside(sigma_sak,w_sak-w_sak_pka*PKA)-Sak),
    # Master regulators
    D(TORC1) ~ gammator*(TORC1_T*soft_Heaviside(sigma_tor,w_torc_glut*Glutamine+w_torc_ego*EGO-w_torc_egoin*(1-EGO)-w_torc-w_torc_snf*Snf1)-TORC1),
    D(Snf1) ~ gammasnf*(Snf1_T*soft_Heaviside(sigma_snf,-w_snf_glc*Carbon+w_snf_sak*Sak-w_snf)-Snf1),
    D(PKA) ~ gammapka*(PKA_T*soft_Heaviside(sigma_pka, w_pka_camp*cAMP-w_pka-w_pka_sch9*Sch9)-PKA),
    D(Sch9) ~ gammasch9*(Sch9_T*soft_Heaviside(sigma_sch9,w_sch9_torc*TORC1-w_sch9)-Sch9),
    # Downstream responses
    D(Gcn2) ~ gamma_gcn2*(Gcn2_T*soft_Heaviside(sigma_gcn2,w_gcn-w_gcn_torc*Sch9)-Gcn2),
    D(Gcn4) ~ (Gcn4_T*soft_Heaviside(sigma_gcn4,w_gcn4_gcn2_trna*minimum([Gcn2, tRNA_sensitivity*(tRNA_total-minimum([tRNA_total, Glutamine]))])-w_gcn4)-Gcn4), # gamma_gcn4 ska läggas in, saknar värde
    D(eIF) ~ gammaeif*(eIF_T*soft_Heaviside(sigma_eif,w_eif-w_eif_gcn2*Gcn2)-eIF),
    D(Gln3) ~ gammagln3*(Gln3_T*soft_Heaviside(sigma_gln,-w_gln3+w_gln_snf*Snf1+ w_gln_sit*(1-TORC1))-Gln3),
    D(Gln1) ~ gammagln1*(Gln1_T*soft_Heaviside(sigma_gln1,w_gln1_gln3*Gln3-w_gln1)-Gln1),
    D(Rtg13) ~ (Rtg13_T*soft_Heaviside(sigma_rtg,-w_rtg_torc*TORC1+w_rtg)-Rtg13), # gamma_rtg13 ska läggas in, saknar värde
    D(Gis1) ~ (Gis1_T*soft_Heaviside(sigma_gis1,-w_gis_pka*PKA-w_gis_sch*Sch9+w_gis)-Gis1), # gamma_gis1 ska läggas in, saknar värde
    D(Mig1) ~ gamma_mig*(Mig1_T*soft_Heaviside(sigma_mig1,w_mig_pka*PKA-w_mig_snf*Snf1+w_mig)-Mig1),
    D(Dot6) ~ (Dot6_T*soft_Heaviside(sigma_dot,-w_dot_sch_pka*Sch9*PKA+w_dot)-Dot6), # gamma_dot6 ska läggas in, saknar värde
    D(Tps1) ~ gammatps*(Tps1_T*soft_Heaviside(sigma_tps, w_tps_pka*(PKA_T-PKA)-w_tps)-Tps1), # 22, Tps1
    D(Trehalase) ~ gammatre*(Trehalase_T*soft_Heaviside(sigma_trehalase, w_tre_pka*PKA-w_tre)-Trehalase), # 23, Trehalase
    D(Protein) ~ k_pr*ATP*minimum([minimum([Rib, eIF]) minimum([tRNA_total, Glutamine])])*Protein,  # 24, Protein
    D(Rib) ~ k_transcription*(1-Dot6)-k_mRNA_degr*Rib, # 25, Rib
    # Model inputs which represent different nutrient conditions values are diefined as initial conditions
    D(Carbon) ~ 0,
    D(ATP) ~ 0,
    D(Glutamine_ext) ~ 0] 

    @named ODE_sys = ODESystem(eqs)

    u_lookup_table = ["Glutamine(t)", "Cyr1(t)", "Ras(t)", "EGO(t)", "EGOGAP(t)", "cAMP(t)",
    "PDE(t)", "Sak(t)", "TORC1(t)", "Snf1(t)", "PKA(t)", "Sch9(t)",
    "Gcn2(t)", "Gcn4(t)", "eIF(t)", "Gln3(t)", "Gln1(t)", "Rtg13(t)",
    "Gis1(t)", "Mig1(t)", "Dot6(t)", "Tps1(t)", "Trehalase(t)",
    "Protein(t)", "Rib(t)", "Carbon(t)", "ATP(t)", "Glutamine_ext(t)"]