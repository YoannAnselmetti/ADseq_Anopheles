# Cedric Chauve, SFU, May 2015
#
# Generates the correspondance species_name species_id from a given set of gene trees

import sys

tree_file = open(sys.argv[1],"r").readlines()
genomes_id={}
species_list=[]
for l in tree_file:
    if l[0]=="(":
        l1=l.rstrip().replace(","," ").replace("("," ").replace(")"," ").replace("["," ").replace("]"," ").split()
        for b in l1:
            if "|" in b:
                gene_name=b.split("|")[0]
                species_name=b.split("|")[1].split(":")[0]
            elif "NHX" in b and not species_name in species_list:
                species_id1=b.split(":")[3]
                if species_id1[0:2]=="S=" and gene_name!="Loss":
                    species_id=species_id1.replace("S=","")
                    genomes_id[species_name]=species_id
                    species_list.append(species_name)

genomes_info_file = open(sys.argv[2],"r").readlines()
for l in genomes_info_file:
    
    l1=l.split()
    species_name=l1[0]
    if species_name in species_list:
        print genomes_id[species_name]+"\t"+l,
