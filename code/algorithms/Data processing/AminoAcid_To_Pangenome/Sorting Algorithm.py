#This function reads a .faa (A gene transcipted in a list of aminoacids)
#It returns a list with the name of the gene and a list with its corresponding amioacid-chain
global correlation_factor
correlation_factor = 0.95
import Smith_Waterman_Algorithm as SWA
import math
import pandas as pd

def read_input(name):
    name_list=[]
    function_list = []
    gene_list = []
    size_list = []
    filedirec = "./PROKKA_" + name + ".faa"
    with open(filedirec,'r') as f:
        #In .faa file each gene starts with the notation ">"
        genome_input = f.read().split(">")
        genome_input.pop(0)
        #Lists cointaning each gene identifier(name), function,  size, and AA chain is created
        for row in range(len(genome_input)):
            gene = genome_input[row].split("\n")
            name_function = gene[0].split(" ")
            name_list.append(name_function[0])
            function_jointed = ' '.join(name_function[1:len(name_function)])
            function_list.append(function_jointed)
            gene = gene[1:len(gene)]
            jointed_gene = ''.join(gene)
            size = len(jointed_gene)
            gene_list.append(jointed_gene)
            size_list.append(size)
            #Make the code work for 100 genes REMOVE TO MAKE CODE WORK FOR ALL THE GENES (not advised-takes a long time)
            if row == 100:
                break

    return gene_list, name_list,function_list, size_list

#This function Verifies if two different genes inputs are similar enougth to be called the same gene
def Verify_Gene_Resemblance(Gene, New_gene,size,newsize):
    #This if statement is used for code optimization (quicker to check size than to create the matrix)
    if size*correlation_factor < newsize <(size*(2-correlation_factor)):
        H, Score, Normalized_Score = SWA.matrix(Gene,New_gene)
        if Normalized_Score > correlation_factor:
            return True
    else:
        return False
        

#This function adds the input gene information to the database
def UpdateDatabase(Pangenome, Gene, size, size_list, identifier, identifier_list, gene_function, function_list):
    Match = False
    if len(Pangenome) > len(identifier_list):
                identifier_list[:] = identifier_list + [math.nan] * (len(Pangenome) - len(identifier_list))
    #Small genes should not be added to the pangenome
    if size > 40:
        for i in range(len(Pangenome)):
            if Verify_Gene_Resemblance(Pangenome[i],Gene,size,size_list[i]):
                identifier_list[i]= identifier
                Match = True
        if not Match:
            Pangenome.append(Gene)
            identifier_list.append(identifier)
            size_list.append(size)
            function_list.append(gene_function)
    return Pangenome, identifier_list, size_list, function_list

#From the Lists made the function outputs a .csv file            
def WritePandasDF_Data(Pangenome, Identifiers, Size_list, Function_list, Bacteria_identifiers):
    for i in range(len(Identifiers)):
        Identifiers[i][:] = Identifiers[i] + [math.nan] * (len(Pangenome) - len(Identifiers[i]))
    Data = [Pangenome] + [Function_list] + [Size_list]
    for i in range(len(Identifiers)):
        Data.append(Identifiers[i])
    df = pd.DataFrame(Data).transpose()
    names = ["Gene","Annotation","Number of AA"] + Bacteria_identifiers
    df = pd.DataFrame({names[i]: col for i, col in enumerate(Data)})
    df.to_csv("Pan_Genome.csv")



def main():
    pan_genome = []
    Identifiers = []
    Size_list = []
    function_list = []
    #Studied bacteria
    bacteria_identifiers = ["10232024","GCF_000007445"]
    #Get each list for each bacteria
    for bacteria in range(len(bacteria_identifiers)):
        new_genome, new_gene_name, new_function_list, new_size_list = read_input(bacteria_identifiers[bacteria])
        identifier_list = []
        #Test each gene and add information to database
        for gene in range(len(new_genome)):
            New_gene = new_genome[gene]
            size = new_size_list[gene]
            identifier = new_gene_name[gene]
            gene_func = new_function_list[gene]
            pan_genome, identifier_list, Size_list, function_list  = UpdateDatabase(pan_genome, New_gene, size, Size_list, identifier, identifier_list, gene_func, function_list)
        #This list of lists cointains information about gene position for each bacteria
        Identifiers.append(identifier_list)
    #Output Data
    WritePandasDF_Data(pan_genome,Identifiers,Size_list,function_list,bacteria_identifiers)

main()


#WritePandasDF_Data(["1","2","3","4"],[[1,2,3,4],[5,6,7]],[1,2,3,4],["a","b","c","d"],["Bact_1","Bact_2"])












