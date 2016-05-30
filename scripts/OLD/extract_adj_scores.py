# Extracting the scores of selected adjacencies into a sorted list

def extract_adj_scores(input_file):

    adj_file=open(input_file,"r").readlines()
    output_list=[]
    
    (nb_adj_read,nb_adj_written)=(0,0)

    for adj in adj_file:
        if adj[0]==">": # Reading a scaffolding adjacency
            adj1=adj.rstrip().split(" ")[1].split(":")
            (c_species,c_gene1,c_gene2)=(adj1[0],adj1[2],adj1[4]) 
            # Species and Genes used to define the contigs adjacency
            c_adj=adj.rstrip()[2:]
            nb_adj_read+=1
        elif adj[0]=="\t": # Looking for the score of the adjacency
            adj1=adj.rstrip().split("\t")
            (species,gene1,gene2,score)=(adj1[1].split(":")[1],adj1[2].split(":")[2],adj1[3].split(":")[2],adj1[4])
            if (species,gene1,gene2)==(c_species,c_gene1,c_gene2) or (species,gene2,gene1)==(c_species,c_gene1,c_gene2):
                output_list.append((float(score),c_adj))
                nb_adj_written+=1

    # Checking that indeed every read adjacency was output with its score
    assert (nb_adj_read==nb_adj_written),"Error: missing written adjacency"

    return(sorted(output_list))
