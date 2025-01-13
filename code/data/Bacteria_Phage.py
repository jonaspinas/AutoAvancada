import pandas as pd
#This file will Create a .csv file containing bacteria phage number

data = {
#Bacteria list
"bacteria_list" : ["GCF_000005845","GCF_000007445","GCF_000013305","GCF_000017985","GCF_000019425","GCF_000025745","GCF_000026245",
            "GCF_000026265","GCF_000026345","GCF_000332755","GCF_000468515","GCF_000599665","GCF_001442495","GCF_001900295","GCF_001900315",
            "GCF_002058765", "GCF_003018255", "GCF_003095635", "GCF_003722195", "GCF_004377995","GCF_008462425","GCF_010092965",
            "GCF_013201505","GCF_902810335"],

#Data Obtained through PHASTEST
"fagos" : [8,8,2,8,9,14,6,3,14,10,5,7,4,5,6,9,9,4,19,4,13,7,1,9]
}


results_df = pd.DataFrame(data)

results_df.to_csv("Phages.csv")
