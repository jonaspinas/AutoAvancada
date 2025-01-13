import pandas as pd
#This file will utilize the gene_function_list obtained in Function_Class_Counter to quantify bacteria genes
# Load the file
url = "https://raw.githubusercontent.com/jonaspinas/Assignment/refs/heads/main/gene_presence_absence.csv"
data = pd.read_csv(url)

#List of gene functions classsification
gene_functions_list = [["hypothetical","protein"],["drug"],["hth-type"],["lipoprotein"],["membrane","protein"],["ribosomal","protein"],["cell division"],["repressor"],
                       ["toxin-antitoxin"],["antitoxin"],["toxin"],["adhesin"],["hemolysin"],["flagellar"],["protease"],
                       ["fimbri"],["pilu"],["host","factor"],["virulence"],["isozyme"],["carrier","protein"],["binding","protein"],
                       ["transport"],["hydrolase"],["kinase"],["oxidase"],["oxidoreductase"],["regulatory","protein"],
                       ["pts system"], ["elongation factor"],["hydroxylase"],["reductase"],["transposase"],["transferase"],
                       ["phosphodiesterase"],["phosphatase"],["ligase"],["aminase"],["synthase"],["lyase"],["phosphorylase"],["helicase"],["crossover"],
                       ["ribonuclease"],["translocase"],["oxygenase"],["polymerase"],["crispr"],["aldolase"],["enzyme"],["mutase"],["putative","protein"]]

#Bacteria list
bacteria_list = ["GCF_000005845","GCF_000007445","GCF_000013305","GCF_000017985","GCF_000019425","GCF_000025745","GCF_000026245",
            "GCF_000026265","GCF_000026345","GCF_000332755","GCF_000468515","GCF_000599665","GCF_001442495","GCF_001900295","GCF_001900315",
            "GCF_002058765", "GCF_003018255", "GCF_003095635", "GCF_003722195", "GCF_004377995","GCF_008462425","GCF_010092965",
            "GCF_013201505","GCF_902810335"]

#Create new df
new_df = []

# Calculate percentage for each word in `most_used_words_list`
for bacteria in bacteria_list:
    gene_total_count = [0] * len(gene_functions_list)
    #Filter for the gene of each bacteria 
    filtered_df = data[data[bacteria].notnull()]
    #Get the genes function
    annotations = filtered_df['Annotation']

    for annotation in annotations:
        annotation = annotation.lower()
        for i in range(len(gene_functions_list)):
            words_counted = 0
            for word in gene_functions_list[i]:
                if annotation.count(word) > 0:
                    words_counted = words_counted + 1
            if words_counted >= len(gene_functions_list[i]):
                gene_total_count[i] = gene_total_count[i] + 1
                break

    new_df.append([bacteria] + gene_total_count)

#Join function words into a single string
for i in range(len(gene_functions_list)):
    gene_functions_list[i] = " ".join(gene_functions_list[i]).capitalize()

results_df = pd.DataFrame(new_df, columns=[["Bacteria"] + gene_functions_list])
results_df.to_csv("Bacteria_Gene_Quantified.csv")
