__author__ = "Cedric Chauve"
__date__ = "May 2016"

# Importing important files from genes data, DeClone and BESST computations.
# This functions do not perform any processing of the data

import sys
from operator import itemgetter

# --------------- Reading tab separated files -----------------------
def read_tab_file(f): # return a list of splitted lines for a tab separated fine of filename f
    f = open(f,"r").readlines()
    o = []
    for l in f:
        if l[0]!="#":
            o.append(l.rstrip().split("\t"))
    return o

# --------------- Reading genes file ---------------------------------
# Genes file format
# Anopheles_albimanus  KB672286  GF0003303  AALB004105  -  461220  466589; TO DO: add exon information

# Genes - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
__GENES_ID=[]     # List of all gene names
__GENES_INFO={}   # Info on each gene
__GENES=[__GENES_ID,__GENES_INFO]

def GENES_list():        # List of all gene names
    return __GENES[0]
def GENE_species(gene):  # species containing gene g given by name
    return __GENES[1][gene][0]
def GENE_scf(gene):      # scaffold/contig/chromosome containing gene g given by name
    return __GENES[1][gene][1]
def GENE_pos(gene):      # realtive position of gene g (by name) in its scaffold/contig/chromosome
    return __GENES[1][gene][6]
def GENE_sign(gene):     # orientation of gene g (by name) in its scaffold/contig/chromosome
    return __GENES[1][gene][3]
def GENE_start(gene):    # start position of  gene g (by name) in its scaffold/contig/chromosome
    return __GENES[1][gene][4]
def GENE_end(gene):      # end position of  gene g (by name) in its scaffold/contig/chromosome
    return __GENES[1][gene][5]
def GENE_og(gene):       # orthogroup of gene g (by name)
    return __GENES[1][gene][2]

# Scaffolds/contigs- - - - - - - - - - - - - - - - - - - - - - - - - - 
__SCF_ID=[]       # ID of scaffolds
__SCF_GENES={}    # Genes on a scaffold
__SCF_SPECIES={}  # Species of a scaffold
__SCF=[__SCF_ID,__SCF_GENES,__SCF_SPECIES]

def SCF_list():         # List of all scaffolds/contigs/chromosomes
    return __SCF[0]
def SCF_genes(scf):     # List of all genes in scaffold/contig/chromosome scf (assumption: no 2 identical names in different species)
    return __SCF[1][scf]
def SCF_nbgenes(scf):   # Number of  genes in scaffold/contig/chromosome scf (assumption: no 2 identical names in different species)
    return len(__SCF[1][scf])
def SCF_species(scf):   # Species containing scaffold/contig/chromosome scf (assumption: no 2 identical names in different species)
    return __SCF[2][scf]

# Translating orientation string +/-/? into integers +1/-1/0
__ORIENT_2_INT={"+":1,"-":-1,"?":0}
__INT_2_ORIENT={1:"+",-1:"-",0:"?"}

# Main import function - - - - - - - - - - - - - - - - - - - - - - - - -
# Read a gene file and populates the structures GENES and SCF
# Sorting genes by species, the chromosome/scaffold/contig, then position
def GENES_import(genes_file): # genes_file is a list corresponding to the splitted lines of a gene file
    genes_file_aux=[]
    for l in genes_file:
        genes_file_aux.append((l[0],l[1],l[2],l[3],l[4],int(l[5]),int(l[6])))
    genes_file_aux.sort(key=itemgetter(5))
    genes_file_aux.sort(key=itemgetter(1))
    genes_file_aux.sort(key=itemgetter(0))
    prev_scf=""
    for (species,scf,og,gene,sign,start,end) in genes_file_aux:
        if scf!=prev_scf:
            __SCF_ID.append(scf)
            __SCF_GENES[scf]=[]            
            position=0
        position+=1
        __GENES[0].append(gene)
        __GENES[1][gene]=(species,scf,og,__ORIENT_2_INT[sign],start,end,position)
        __SCF_GENES[scf].append(gene)                          
        __SCF_SPECIES[scf]=species
        prev_scf=scf        
    __GENES[0].append("NA")
    __GENES[1]["NA"]=("ANCESTRAL","NA","NA",0,0,0,0)
    __SCF_ID.append("NA")
    __SCF_GENES["NA"]=[]
    __SCF_SPECIES["NA"]="NA"

# Orthogroups - - - - - - - - - - - - - - - - - - - - - - - - - - 
__OG_ID=[]      # List of IDs of orthogroups 
__OG_GENES={}   # Orthogroup indexed by genes names
__OG=[__OG_ID,__OG_GENES]

def OG_list():      # List of all orthogroups
    return __OG[0]
def OG_genes(og):   # List of all genes in orthogroup og
    return __OG[1][og]
def OG_size(og):    # Number of genes in orthogroup og
    return len(__OG[1][og])

# Importing the orthogroups
def OG_import(genes_file): # gene_file is a list corresponding to the splitted lines of a gene file
    genes_file_aux=[]
    for l in genes_file:
        genes_file_aux.append((l[0],l[1],l[2],l[3],l[4],int(l[5]),int(l[6])))
    genes_file_aux.sort(key=itemgetter(2))
    prev_og=""
    for (species,scf,og,gene,sign,start,end) in genes_file_aux:
        if og!=prev_og:
            __OG_ID.append(og)
            __OG_GENES[og]=[]            
        __OG_GENES[og].append(gene)                          
        prev_og=og

# Genomes - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
__GENOMES_ID=[]    # Names of genomes
__GENOMES_SCF={}   # Scaffolds of a given genome
__GENOMES=[__GENOMES_ID,__GENOMES_SCF]

def GENOMES_list():     # List of all genomes names
    return __GENOMES[0]
def GENOMES_scf(g):     # List of all scaffols/contigs/chromosomes in a genome given by name
    return __GENOMES[1][g]
def GENOMES_nbscf(g):   # Number of scaffolds/contigs/chromosomes in a genome given by name
    return len(__GENOMES[1][g])
def GENOMES_nbgenes(g): # Number of genes in a genome given by name
    res=0
    for scf in GENOMES_scf(g):
        res+=SCF_nbgenes(scf)
    return res

# Importing genomes
def GENOMES_import(genes_file): # gene_file is a list corresponding to the splitted lines of a genes file
    genes_file_aux=[]
    for l in genes_file:
        genes_file_aux.append((l[0],l[1]))
    genes_file_aux.sort(key=itemgetter(1)) # sorting by scaffold name 
    genes_file_aux.sort(key=itemgetter(0)) # then by species name
    prev_species=""
    prev_scf=""
    for (species,scf) in genes_file_aux:
        if species!=prev_species:
            __GENOMES_ID.append(species)
            __GENOMES_SCF[species]=[]            
        if scf!=prev_scf:
            __GENOMES_SCF[species].append(scf)                          
        prev_scf=scf
        prev_species=species

# --------------- Reading besst file ---------------------------------
# BESST file format
#    Anopheles_albimanus KB672287 KB672445 + + -747.187384203 GF0006744 GF0010639 AALB003904 AALB007558 - + 21633.8126158 0.372384235922 0.362637362637 91

__BESST_LINKS={} # BESST links indexed by pairs of contigs and pairs of genes
__BESST_NGBS={}  # Neighbours of a gene or contig
__BESST=[__BESST_LINKS,__BESST_NGBS]

def BESST_id(pair): # pair = (i1,i2) = pair of genes or pair of contigs/scaffold/chromosome: returns the corresp. besst entry
    return __BESST[0][pair]
def BESST_species(l):  # species of a besst entry
    return l[0]  
def BESST_ctg1(l):     # contig 1 of a besst entry
    return l[1]
def BESST_ctg2(l):     # contig 2 of a besst entry
    return l[2]
def BESST_orctg1(l):   # orientation of ctg 1 in a besst entry
    return l[3]
def BESST_orctg2(l):   # orientation of ctg 2 in a besst entry
    return l[4]
def BESST_distc(l):    # distance between contigs in a besst entry
    return float(l[5])
def BESST_gene1(l):    # gene 1 in a besst entry
    return l[6]
def BESST_gene2(l):    # gene 2 in a besst entry
    return l[7]
def BESST_distg(l):    # distance between genes in a besst entry
    return float(l[8])
def BESST_vscore(l):   # besst vscore
    return float(l[9])
def BESST_dscore(l):   # besst dscore
    return float(l[10])
def BESST_nbpe(l):     # number of paired-end reads supporting a besst entry
    return int(l[11])
def BESST_adj(i1,i2):  # test if i2 is the mate of i1 to define a besst entry: boolean
    return i2 in __BESST[1][i1]

# Importing a besst file
def BESST_import(besst_file): # besst_file is a list of the splitted lines of a BESST file
    for scf in SCF_list():
        __BESST_NGBS[scf]=[]
    for g in GENES_list():
        __BESST_NGBS[g]=[]
    for (s,c1,c2,oc1,oc2,d1,gf1,gf2,g1,g2,og1,og2,d2,vs,ds,nbl) in besst_file:
        __BESST_LINKS[(g1,g2)]=(s,c1,c2,__ORIENT_2_INT[oc1],__ORIENT_2_INT[oc2],d1,g1,g2,d2,vs,ds,nbl)
        __BESST_LINKS[(g2,g1)]=(s,c2,c1,-1*__ORIENT_2_INT[oc2],-1*__ORIENT_2_INT[oc1],d1,g2,g1,d2,vs,ds,nbl)
        __BESST_LINKS[(c1,c2)]=(s,c1,c2,__ORIENT_2_INT[oc1],__ORIENT_2_INT[oc2],d1,g1,g2,d2,vs,ds,nbl)
        __BESST_LINKS[(c2,c1)]=(s,c2,c1,-1*__ORIENT_2_INT[oc2],-1*__ORIENT_2_INT[oc1],d1,g2,g1,d2,vs,ds,nbl)
        __BESST_NGBS[g1].append(g2)
        __BESST_NGBS[g2].append(g1)
        __BESST_NGBS[c1].append(c2)
        __BESST_NGBS[c2].append(c1)

# --------------- Reading DeClone results ---------------------------------------------

__DECLONE_ID=[]           # List of instances IDs
__DECLONE_INSTANCES={}    # List of adjacencies organized by instances
__DECLONE_ADJACENCIES=[]  # List of all adjacencies
__DECLONE=[__DECLONE_ID,__DECLONE_INSTANCES,__DECLONE_ADJACENCIES]

def DECLONE_instances_list():    # list of instances IDs
    return __DECLONE[0]
def DECLONE_adj_per_instance(i): # list of adjacencies in an instance
    return __DECLONE[1][i]
def DECLONE_adjacencies_list():  # list of all adjacencieds
    return __DECLONE[2]
def DECLONE_adj_species_id(adj):    # species ID of an adjacency
    return adj[0][0]
def DECLONE_adj_species_name(adj):  # species name of an adjacency
    return adj[0][1]
def DECLONE_adj_gene1_tree(adj):    # tree ID of the gene 1 of an adjacency
    return adj[1][0]
def DECLONE_adj_gene1_node(adj):    # node ID in its tree of gene 1 of an adjacency
    return adj[1][1]
def DECLONE_adj_gene1_name(adj):    # name of gene 1 in an adjacency (NA if ancestral)
    return adj[1][2]
def DECLONE_adj_gene2_tree(adj):    # tree ID of the gene 2 of an adjacency
    return adj[2][0]
def DECLONE_adj_gene2_node(adj):    # node ID in its tree of gene 2 of an adjacency
    return adj[2][1]
def DECLONE_adj_gene2_name(adj):    # name of gene 2 in an adjacency (NA if ancestral)
    return adj[2][2]
def DECLONE_adj_score(adj):         # DeClone score of an adjacency
    return adj[3]
def DECLONE_adj_instance(adj):      # instance containing an adjacency
    return adj[4]
def DECLONE_adj_besst(adj):         # besst support for an adjacency (BESST format)
    return adj[5]
def DECLONE_adj_orctg1(adj):        # orientation of contig 1 in an adjacency
    return adj[6][0]
def DECLONE_adj_orctg2(adj):        # orientation of contig 2 in an adjacency
    return adj[6][1]

def __split_adjacency(l1): # Local: Splitting an adjacency into (????)
    l11=l1[1].split(":")
    l12=l1[2].split(":")
    l13=float(l1[3])
    l14=l1[4].split(":")
    l15=l1[5].split("|")
    return (l1[0],l11,l12,l13,l14,l15)
                  
# Orienting a scaffold 
# Approach 1: looking at the given extremity: only scaffolds of size 1 can not be oriented
# Approach 2: looking if the observed homologous adjacencies all agree on an orientation
# c is a correction: 1 for first gene of an adjacency, -1 for second gene
def __orient_scf(gene,instance,c): # Local: Orienting a scaffold containing a gene in an instance
    (nbg,pos,sign)=(SCF_nbgenes(GENE_scf(gene)),GENE_pos(gene),GENE_sign(gene))
    result=0
    if sign!=0 and nbg>1 and pos==nbg:
        result=c
    elif sign!= 0 and nbg>1 and pos==1:
        result=-c
    elif sign!=0: # gene is on a scaffold containing one gene
        observed={1:0,-1:0} # number of observed extant adjacencies indication + and - respectively
        for adj in DECLONE_adj_per_instance(instance):
            g1=DECLONE_adj_gene1_name(adj)
            g2=DECLONE_adj_gene2_name(adj)
            (pos1,scf1,sign1)=(GENE_pos(g1),GENE_scf(g1),GENE_sign(g1))
            (pos2,scf2,sign2)=(GENE_pos(g2),GENE_scf(g2),GENE_sign(g2))
            signs={1:sign1,-1:sign2}
            if DECLONE_adj_score(adj)==1 and scf1==scf2 and abs(pos1-pos2)==1:
                # observed adjacency between extant genes
                if pos1<pos2 and signs[c]==1:
                    observed[1*sign]+=1
                elif pos1<pos2 and signs[c]==-1:
                    observed[-1*sign]=+1
                elif pos1>pos2 and signs[c]==1:
                    observed[-1*sign]+=1
                elif pos1>pos2 and signs[c]==-1:
                    observed[1*sign]+=1
        if observed[1]>0 and observed[-1]==0:
            result=1
        elif observed[-1]>0 and observed[1]==0:
            result=-1
    return result

# Importing DeClone adjacencies
def DECLONE_import(declone_file):
    # Step 1: no orientation
    __DECLONE_INSTANCES_AUX={}
    declone_file.sort(key=itemgetter(5))
    prev_instance=("","","")    
    for l in declone_file:
        (l10,l11,l12,l13,l14,l15)=__split_adjacency(l)
        gene1=l14[1]
        gene2=l14[2]
        adj=((l10,l14[0]),(l11[0],l11[1],l14[1]),(l12[0],l12[1],l14[2]),l13,(l15[0],l15[1],l15[2]),("NA","NA","NA","NA","NA")) 
        tl15=tuple(l15)
        if tl15!=prev_instance:
            __DECLONE_ID.append(tl15)
            __DECLONE_INSTANCES_AUX[tl15]=[]
        __DECLONE_INSTANCES_AUX[tl15].append(adj)
        prev_instance=tl15
    # Step 2: Orienting DeClone adjacencies
    for instance in DECLONE_instances_list():
        __DECLONE_INSTANCES[instance]=[]
        for adj in __DECLONE_INSTANCES_AUX[instance]:
            gene1=DECLONE_adj_gene1_name(adj)
            gene2=DECLONE_adj_gene2_name(adj)
            if DECLONE_adj_score(adj)<1:
                or1=__orient_scf(gene1,instance,1)
                or2=__orient_scf(gene2,instance,-1)
            else:
                (or1,or2)=(0,0)
            __DECLONE_INSTANCES[instance].append((adj[0],adj[1],adj[2],adj[3],adj[4],adj[5],(or1,or2)))
            __DECLONE_ADJACENCIES.append((adj[0],adj[1],adj[2],adj[3],adj[4],adj[5],(or1,or2)))
