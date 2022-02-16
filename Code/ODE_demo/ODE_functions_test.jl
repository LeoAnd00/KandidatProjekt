# 1=Glutamin, Cyr1=2, 3=Ras, 4=EGO, 5= EGOGAP, 6= cAMP
#  7=PDE, 8=Sak 9=TORC1, 10=Snf1, 11=PKA, 12=Sch9
# 13=Gcn2, 14=Gcn4, 15=eIF, 16=Gln3, 17=Gln1, 18=Rtg13
# 19=Gis1, 20=Mig1, 21=Dot6 22=Tps1, 23=Trehalase
# 24=Protein, 25=Rib

include("Get_parameter_values_test.jl")

function soft_Heaviside(sigma, W)
    1/(1+exp(-sigma*W))
end

function ODE_functions(du, u, t)
    du[22] = Tps1_parameters["gammatps"](*u[22]*soft_Heaviside(Tps1_parameters["sigma_tps"], Tps1_parameters["w_tps_pka"]*(u[?] - u[11])))
    du[23] = Trehalase_parameters["gammatre"]*(u[?]*soft_Heaviside(Trehalase_parameters["sigma_trehalase"], Trehalase_parameters["w_tre_pka"]*u[11]-Trehalase_parameters["w_tre"])-u[23])
    du[24] = Protein_parameters["k_pr"]*Protein_parameters["ATP"]*minimum([minimum([u[25] u[15]]) minimum(Protein_parameters["tRNA_total"] u[1])])*u[24]
    du[25] = Rib_parameters["k_transcription"]*(1-u[21])-Rib_parameters["k_mRNA_degr"]*u[25]
end

