# Cedric Chauve, SFU, May 2015
#
# From a pair of files Trees+adjacencies, generates a list of pairs trees+adjacencies files and an instances list file
import sys 

# argv[1] adjacenciy classes with gene trees
# argv[2] extant adjacencies
# argv[3] output_directory
# argv[4] list of instances
# argv[4]+"_relevant" list of relevant instances for extant scaffolding


output_directory=sys.argv[3]
list_file=open(sys.argv[4],"w")
list_file_relevant=open(sys.argv[4]+"_relevant","w")

instances_list = []
instances_genes = {}  # list of genes for each instance (gene number, species number, tree id)
adjacencies = {}      # list of adjacencies
degrees = {}          # degree of each gene                                                              
species_list = []     # list of species
                  
print "READING AND WRITING TREE FILES"
trees_file = open(sys.argv[1],"r").xreadlines()
for l in trees_file:
    if l[0]=="A": # new adjacency class
        instance_id=l.rstrip().split()[2]  # label/ID of the instance
        instances_list.append(instance_id)
        instance_trees_file=open(output_directory+"/"+instance_id.replace("|","_")+"_trees","w") # creating the corresponding instance file
        list_file.write(l)
        instances_genes[instance_id]=[]
        tree_id=0
    elif l[0]=="G":
        tree_id+=1
    elif l[0]=="(":
        t=l.rstrip().replace("("," ").replace(":"," ").replace(","," ").split()
        for t1 in t:
            if "|" in t1 and t1[0:4]!="Loss":
                gene=t1.split("|")[0]
                species=t1.split("|")[1]                
                if species not in species_list:
                    species_list.append(species)
                adjacencies[gene] = []
                degrees[gene]=0 # ancestral gene
                instances_genes[instance_id].append((gene,species,tree_id))
    instance_trees_file.write(l)
    #if l[0]=="*":
    #    instance_trees_file.close()

print "READING ADJACENCIES FILE"
adjacencies_file = open(sys.argv[2],"r").readlines()
for l in adjacencies_file:
    l1=l.rstrip().split()
    adjacencies[l1[0]]=[]
    adjacencies[l1[1]]=[]
    degrees[l1[0]]=0
    degrees[l1[1]]=0

for l in adjacencies_file:
    l1=l.rstrip().split()
    adjacencies[l1[0]].append(l1[1])
    adjacencies[l1[1]].append(l1[0])
    degrees[l1[0]]+=1
    degrees[l1[1]]+=1

print "WRITING ADJACENCIES FILES"
degree_above={}
degree_above[1]={}
degree_above[2]={}
for instance_id in instances_list:
    instance_adj_file=open(output_directory+"/"+instance_id.replace("|","_")+"_adj","w")
    genes_list=instances_genes[instance_id]
    for s in species_list:                                                         
        degree_above[1][s]=0
        degree_above[2][s]=0
    for gl in genes_list:
        g1=gl[0]
        for g2 in adjacencies[g1]:
            if (g2 not in genes_list or g2>g1):
                instance_adj_file.write(g1+" "+g2+"\n")
        if degrees[gl[0]]<2:
            degree_above[gl[2]][gl[1]]+=1    
    relevant=False
    for s in species_list: 
        if degree_above[1][s]>=1 and degree_above[2][s]>=1:
            relevant=True
    if relevant:
        list_file_relevant.write("Adjacency class "+instance_id+"\n")
    instance_adj_file.close()

    

#Adjacency class 0|0|1082
#
#Gene tree 0|1048
#(((((((((ENSGGOP00000014946|Gorilla_gorilla:1[&&NHX:Ev=Extant:ND=954:S=20],Loss|19:1[&&NHX:Ev=GLos:ND=955:S=19]):1[&&NHX:Ev=Spec:ND=956:S=21],Loss|Pongo_abelii:1[&&NHX:Ev=GLos:ND=957:S=22]):1[&&NHX:Ev=Spec:ND=958:S=23],((ENSPANP00000004733|Papio_anubis:1[&&NHX:Ev=Extant:ND=959:S=12],ENSMMUP00000037185|Macaca_mulatta:1[&&NHX:Ev=Extant:ND=960:S=13]):1[&&NHX:Ev=Spec:ND=961:S=14]
