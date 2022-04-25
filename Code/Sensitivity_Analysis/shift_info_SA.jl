include("../../Data/exp_data.jl")

# 1=Glutamine, Cyr1=2, 3=Ras, 4=EGO, 5=EGOGAP, 6=cAMP
# 7=PDE, 8=Sak 9=TORC1, 10=Snf1, 11=PKA, 12=Sch9
# 13=Gcn2, 14=Gcn4, 15=eIF, 16=Gln3, 17=Gln1, 18=Rtg13
# 19=Gis1, 20=Mig1, 21=Dot6 22=Tps1, 23=Trehalase
# 24=Protein, 25=Rib

out_indices = [10,6,6,12,12,12,12,12,20,25]

all_tvals = [t_Snf1, t_cAMP, t_sch9Delta_cAMP, t_Sch9P_glutamine_H, t_Sch9P_glutamine_L,
t_Sch9_gtr1Delta, t_Sch9_glucose_starve, t_Sch9_glucose_relief, t_Mig1_glucose_relief, t_Rib_rap]

pre_shifts = [(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 0.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)]

post_shifts = [(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 0.3),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 0.0, ATP => 0.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0),
(Carbon => 1.0, ATP => 1.0, Glutamine_ext => 1.0)]

all_pre_p_const_changes = [[],[],Dict("Sch9_T" => (Sch9_T => 0.0)),[],[],
Dict("EGO_T" => (EGO_T => 0.0)),[],[],[],[]]

all_post_p_const_changes = [[],[],Dict("Sch9_T" => (Sch9_T => 0.0)),[],[],[],[],[],[],
Dict("TORC1_T" => (TORC1_T => 0.0))]

all_pre_p_var_changes = [[],[],[],[],[],Dict("w_torc_ego" => 0.0,
"w_torc_egoin" => 0.0),[],[],[],[]]