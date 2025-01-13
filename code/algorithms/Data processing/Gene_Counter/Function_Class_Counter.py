from collections import Counter
import pandas as pd

#This code was used to create the different function classification, in order to quantify gene presence in bacteria.
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
# Load the file
url = "https://raw.githubusercontent.com/jonaspinas/Assignment/refs/heads/main/gene_presence_absence.csv"
data = pd.read_csv(url)

# Extract the 'Annotation' column
annotations = data['Annotation']
nested_word_vector = []
gene_functions_list = [["hypothetical","protein"],["drug"],["hth-type"],["lipoprotein"],["membrane","protein"],["ribosomal","protein"],["cell division"],["repressor"],
                       ["toxin-antitoxin"],["antitoxin"],["toxin"],["adhesin"],["hemolysin "],["flagellar"],["hemolysin"],["protease"],
                       ["fimbri"],["pilu"],["host","factor"],["virulence"],["isozyme"],["carrier","protein"],["binding","protein"],
                       ["transport"],["hydrolase"],["kinase"],["oxidase"],["oxidoreductase"],["regulatory","protein"],
                       ["pts system"], ["elongation factor"],["hydroxylase"],["reductase"],["transposase"],["transferase"],
                       ["phosphodiesterase"],["phosphatase"],["ligase"],["aminase"],["synthase"],["lyase"],["phosphorylase"],["helicase"],["crossover"],
                       ["ribonuclease"],["translocase"],["oxygenase"],["polymerase"],["crispr"],["aldolase"],["enzyme"],["mutase"],["putative","protein"]]
gene_total_count = [0] * len(gene_functions_list)

# Split each annotation into words and verify if the annotation matches a function
for annotation in annotations:
    annotation = annotation.lower()
    Match = False
    for i in range(len(gene_functions_list)):
        words_counted = 0
        for word in gene_functions_list[i]:
            if annotation.count(word) > 0:
                words_counted = words_counted + 1
        if words_counted >= len(gene_functions_list[i]):
            Match = True
            gene_total_count[i] = gene_total_count[i] + 1
            break
    #If a gene doesn't belong to a category, its words will be added to nested_word_vector
    #This list is usefull to create new function classification based on not classified genes
    if not Match:
        name_list = annotation.split() 
        # Filter for words longer then 1 character
        valid_names = [word for word in name_list if len(word) > 1]
        nested_word_vector.append(valid_names)

# Count occurrences of each word globally
flat_word_list = [word for sublist in nested_word_vector for word in sublist]
global_word_counts = Counter(flat_word_list)

# Filter words that appear more than once globally
word_counts_fil = {word: count for word, count in global_word_counts.items() if count > 1}

# Get unique words to analyze
most_used_words_list = list(word_counts_fil.keys())

# Display the results
for i in range(len(gene_total_count)):
    print(gene_functions_list[i],gene_total_count[i])

print(sum(gene_total_count))