# Labeling the genetres according to Anopheles gambiae chromosomes
# X    = all genes on X chromosome
# noX  = no gene on X chromosome
# ambX = genes on X and autosome

import sys
from Bio import Phylo

gene_file=open(sys.argv[1],"r").readlines()
gtrees_file=sys.argv[2]
output_file=open(sys.argv[3],"w")

genes_chr={}

for l in gene_file:
    l1=l.split("\t")
    species=l1[0]
    chrom=l1[1]
    og=l1[2]
    gene=l1[3]
    if species=="Anopheles_gambiae" and gene!="Loss":
        genes_chr[gene]=chrom
    

gtrees=Phylo.parse(gtrees_file,"newick")
status={} # 0 = no info, 1 = X genes, 2 = X+losses, 3 = autosomes, 4 = autosomes+loss, 5 = ambiguous
tree=0 # Tree ID
for gtree in gtrees:
    tree_ID=str(tree)
    clades=gtree.find_clades()
    nbXs=0
    nbAutosomes=0
    nbLosses=0
    for c in clades:
        if c.name!=None and "Anopheles_gambiae" in c.name:
            gene=c.name.split("|")[0]
            if gene=="Loss":
                nbLosses+=1
            elif genes_chr[gene]=="X":
                nbXs+=1
            else:
                nbAutosomes+=1
    if nbXs==0:
        if nbLosses==0:
            status[tree_ID]=3
        else:
            status[tree_ID]=4
    else: #nbXs>0
        if nbAutosomes==0:
            if nbLosses==0:
                status[tree_ID]=1
            else:
                status[tree_ID]=2
        else:
            status[tree_ID]=5        
    tree+=1


status2str={0:"unknown", 1:"chrom.X", 2:"x+losses", 3:"autosome", 4:"aut+losses", 5:"ambiguous"}
output_file.write("# Status of anopheles gambia gene trees regarding their location: chrom.X, autosome, losses, ambiguous\n")
stats={1:0,2:0,3:0,4:0,5:0}
for tree in status.keys():
    if status[tree]>0:        
        output_file.write(tree+"\t"+status2str[status[tree]]+"\n")
        stats[status[tree]]+=1
print "chrom.X\t"+str(stats[1])+"\tchrX+losses\t"+str(stats[2])+"\tautosomes\t"+str(stats[3])+"\taut+losses\t"+str(stats[4])+"\tambiguous\t"+str(stats[5])
