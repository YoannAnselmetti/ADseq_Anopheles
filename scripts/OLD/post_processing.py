# Postprocessing the DeClone/BESST output.
# Input:
# file 1 = genes information (from GFF format + orthogroups)
# file 2 = BESST results
# file 3 = DeClone results

import sys
from operator import itemgetter

# --------------- Reading genes file ---------------------------------
# Genes file format
# Anopheles_albimanus  KB672286  GF0003303  AALB004105  -  461220  466589; TO DO: add exon information

genes_file=[]

GENES_ID=[]     # List of all gene names
GENES_INFO={}   # Info on each gene
GENES=[GENES_ID,GENES_INFO]

def GENES_list():
    return GENES[0]
def GENE_species(gene):
    return GENES[1][gene][0]
def GENE_scf(gene):
    return GENES[1][gene][1]
def GENE_pos(gene):
    return GENES[1][gene][6]
def GENE_sign(gene):
    return GENES[1][gene][3]
def GENE_start(gene):
    return GENES[1][gene][4]
def GENE_end(gene):
    return GENES[1][gene][5]
def GENE_og(gene):
    return GENES[1][gene][2]

SCF_ID=[]       # ID of scaffolds
SCF_GENES={}    # Genes on a scaffold
SCF_SPECIES={}  # Species of a scaffold
SCF=[SCF_ID,SCF_GENES,SCF_SPECIES]

def SCF_list():
    return SCF[0]
def SCF_genes(scf):
    return SCF[1][scf]
def SCF_nbgenes(scf):
    return len(SCF[1][scf])
def SCF_species(scf):
    return SCF[2][scf]

ORIENT_2_INT={"+":1,"-":-1,"?":0}
INT_2_ORIENT={1:"+",-1:"-",0:"?"}

def GENES_import():
    # Sorting genes by (1) species, (2) chromosome/scaffold/contig, (6) position
    genes_file.sort(key=itemgetter(5))
    genes_file.sort(key=itemgetter(1))
    genes_file.sort(key=itemgetter(0))
    prev_scf=""
    for (species,scf,og,gene,sign,start,end) in genes_file:
        if scf!=prev_scf:
            SCF_ID.append(scf)
            SCF_GENES[scf]=[]            
            position=0
        position+=1
        GENES[0].append(gene)
        GENES[1][gene]=(species,scf,og,ORIENT_2_INT[sign],start,end,position)
        SCF_GENES[scf].append(gene)                          
        SCF_SPECIES[scf]=species
        prev_scf=scf        
    GENES[0].append("NA")
    GENES[1]["NA"]=("ANCESTRAL","NA","NA",0,0,0,0)
    SCF_ID.append("NA")
    SCF_GENES["NA"]=[]
    SCF_SPECIES["NA"]="NA"

def GENE_write(g11):
    return g11+":"+GENE_scf(g11)+"_"+str(SCF_nbgenes(GENE_scf(g11)))+":"+str(GENE_pos(g11))+":"+str(GENE_start(g11))+":"+str(GENE_end(g11))+":"+INT_2_ORIENT[GENE_sign(g11)]


OG_ID=[]      # List of IDs of orthogroups 
OG_GENES={}   # Orthogroup indexed by genes names
OG=[OG_ID,OG_GENES]

def OG_list():
    return OG[0]
def OG_genes(og):
    return OG[1][og]
def OG_size(og):
    return len(OG[1][og])

def OG_import():
    genes_file.sort(key=itemgetter(2))
    prev_og=""
    for (species,scf,og,gene,sign,start,end) in genes_file:
        if og!=prev_og:
            OG_ID.append(og)
            OG_GENES[og]=[]            
        OG_GENES[og].append(gene)                          
        prev_og=og

GENOMES_ID=[]    # Names of genomes
GENOMES_SCF={}   # Scaffolds of a given genome
GENOMES=[GENOMES_ID,GENOMES_SCF]

def GENOMES_list():
    return GENOMES[0]
def GENOMES_scf(g):
    return GENOMES[1][g]
def GENOMES_nbscf(g):
    return len(GENOMES[1][g])
def GENOMES_nbgenes(g):
    res=0
    for scf in GENOMES_scf(g):
        res+=SCF_nbgenes(scf)
    return res

def GENOMES_import():
    genes_file.sort(key=itemgetter(1))
    genes_file.sort(key=itemgetter(0))
    prev_species=""
    prev_scf=""
    for (species,scf,og,gene,sign,start,end) in genes_file:
        if species!=prev_species:
            GENOMES_ID.append(species)
            GENOMES_SCF[species]=[]            
        if scf!=prev_scf:
            GENOMES_SCF[species].append(scf)                          
        prev_scf=scf
        prev_species=species

# --------------- Reading besst file ---------------------------------
# BESST file format
#    Anopheles_albimanus KB672287 KB672445 + + -747.187384203 GF0006744 GF0010639 AALB003904 AALB007558 - + 21633.8126158 0.372384235922 0.362637362637 91

besst_file=[]

BESST_LINKS={} # BESST links indexed by pairs of contigs and pairs of genes
BESST_NGBS={}
BESST=[BESST_LINKS,BESST_NGBS]

def BESST_id((i1,i2)):
    return BESST[0][(i1,i2)]
def BESST_species(l):
    return l[0]
def BESST_ctg1(l):
    return l[1]
def BESST_ctg2(l):
    return l[2]
def BESST_orctg1(l):
    return l[3]
def BESST_orctg2(l):
    return l[4]
def BESST_distc(l):
    return float(l[5])
def BESST_gene1(l):
    return l[6]
def BESST_gene2(l):
    return l[7]
def BESST_distg(l):
    return float(l[8])
def BESST_vscore(l):
    return float(l[9])
def BESST_dscore(l):
    return float(l[10])
def BESST_nbpe(l):
    return int(l[11])
def BESST_adj(i1,i2):
    return i2 in BESST[1][i1]

def BESST_import():
    for scf in SCF_list():
        BESST_NGBS[scf]=[]
    for g in GENES_list():
        BESST_NGBS[g]=[]
    for (s,c1,c2,oc1,oc2,d1,gf1,gf2,g1,g2,og1,og2,d2,vs,ds,nbl) in besst_file:
        BESST_LINKS[(g1,g2)]=(s,c1,c2,ORIENT_2_INT[oc1],ORIENT_2_INT[oc2],d1,g1,g2,d2,vs,ds,nbl)
        #BESST_LINKS[(g2,g1)]=(s,c1,c2,ORIENT_2_INT[oc1],ORIENT_2_INT[oc2],d1,g1,g2,d2,vs,ds,nbl)
        BESST_LINKS[(g2,g1)]=(s,c2,c1,-1*ORIENT_2_INT[oc2],-1*ORIENT_2_INT[oc1],d1,g2,g1,d2,vs,ds,nbl)
        BESST_LINKS[(c1,c2)]=(s,c1,c2,ORIENT_2_INT[oc1],ORIENT_2_INT[oc2],d1,g1,g2,d2,vs,ds,nbl)
        #BESST_LINKS[(c2,c1)]=(s,c1,c2,ORIENT_2_INT[oc1],ORIENT_2_INT[oc2],d1,g1,g2,d2,vs,ds,nbl)        
        BESST_LINKS[(c2,c1)]=(s,c2,c1,-1*ORIENT_2_INT[oc2],-1*ORIENT_2_INT[oc1],d1,g2,g1,d2,vs,ds,nbl)        
        BESST_NGBS[g1].append(g2)
        BESST_NGBS[g2].append(g1)
        BESST_NGBS[c1].append(c2)
        BESST_NGBS[c2].append(c1)

def BESST_write(l):
    return BESST_ctg1(l)+":"+INT_2_ORIENT[BESST_orctg1(l)]+":"+BESST_ctg2(l)+":"+INT_2_ORIENT[BESST_orctg2(l)]+":"+str(BESST_distc(l))+":"+str(BESST_distg(l))+":"+str(BESST_vscore(l))+":"+str(BESST_dscore(l))+":"+str(BESST_nbpe(l))

# --------------- Reading DeClone results ---------------------------------------------

declone_file=[]

DECLONE_ID=[]           # List of instances IDs
DECLONE_INSTANCES={}    # List of all instances
DECLONE=[DECLONE_ID,DECLONE_INSTANCES]#,DECLONE_GENE_PAIRS,DECLONE_SPECIES,DECLONE_SCF]

def DECLONE_instances():
    return DECLONE[0]
def DECLONE_instance_write(i):
    return i[0]+":"+i[1]+":"+i[2]

def DECLONE_adj_per_instance(i):
    return DECLONE[1][i]

def DECLONE_adj_species_id(adj):
    return adj[0][0]
def DECLONE_adj_species_name(adj):
    return adj[0][1]
def DECLONE_adj_gene1_tree(adj):
    return adj[1][0]
def DECLONE_adj_gene1_node(adj):
    return adj[1][1]
def DECLONE_adj_gene1_name(adj):
    return adj[1][2]
def DECLONE_adj_gene2_tree(adj):
    return adj[2][0]
def DECLONE_adj_gene2_node(adj):
    return adj[2][1]
def DECLONE_adj_gene2_name(adj):
    return adj[2][2]
def DECLONE_adj_score(adj):
    return adj[3]
def DECLONE_adj_instance(adj):
    return adj[4]
def DECLONE_adj_besst(adj):
    return adj[5]
def DECLONE_adj_orctg1(adj):
    return adj[6][0]
def DECLONE_adj_orctg2(adj):
    return adj[6][1]

def split_adjacency(l1):
    l11=l1[1].split(":")
    l12=l1[2].split(":")
    l13=float(l1[3])
    l14=l1[4].split(":")
    l15=l1[5].split("|")
    return (l1[0],l11,l12,l13,l14,l15)

def adjacency_write(adj):
    species1=DECLONE_adj_species_id(adj)
    species2=DECLONE_adj_species_name(adj)
    species=species1+":"+species2
    gene1=DECLONE_adj_gene1_name(adj)
    gene2=DECLONE_adj_gene2_name(adj)
    (scf1,scf2)=(GENE_scf(gene1),GENE_scf(gene2))
    (nbscf1,nbscf2)=(SCF_nbgenes(scf1),SCF_nbgenes(scf2))
    (or1,og1,pos1,or2,og2,pos2)=(GENE_sign(gene1),GENE_og(gene1),GENE_pos(gene1),GENE_sign(gene2),GENE_og(gene2),GENE_pos(gene2))
    if BESST_adj(gene1,gene2):
        besst=BESST_id((gene1,gene2))
        strbesst=str(BESST_vscore(besst))+":"+str(BESST_dscore(besst))+":"+str(BESST_nbpe(besst))+":"+INT_2_ORIENT[BESST_orctg1(besst)]+":"+INT_2_ORIENT[BESST_orctg2(besst)]+":"+str(int(BESST_distc(besst)))+":"+str(int(BESST_distg(besst)))
    else:
        strbesst="NA:NA:NA:NA:NA:NA:NA"
    strgene1=scf1+":"+str(nbscf1)+":"+gene1+":"+str(pos1)+":"+og1+":"+INT_2_ORIENT[or1]
    strgene2=scf2+":"+str(nbscf2)+":"+gene2+":"+str(pos2)+":"+og2+":"+INT_2_ORIENT[or2]
    instance=DECLONE_adj_instance(adj)
    strinstance=DECLONE_instance_write(instance)
    return species+"\t"+strgene1+"\t"+strgene2+"\t"+str(DECLONE_adj_score(adj))+"\t"+strinstance+"\t"+strbesst


def adjacency_write_short(adj): # Assumption: extant adjacency
    species=DECLONE_adj_species_name(adj)
    (gene1,gene2)=(DECLONE_adj_gene1_name(adj),DECLONE_adj_gene2_name(adj))
    (scf1,scf2)=(GENE_scf(gene1),GENE_scf(gene2))
    (or1,or2)=(DECLONE_adj_orctg1(adj),DECLONE_adj_orctg2(adj))    
    return species+":"+scf1+":"+gene1+":"+scf2+":"+gene2+":"+INT_2_ORIENT[or1]+":"+INT_2_ORIENT[or2]

# Orienting a scaffold 
# Approach 1: looking at the given extremity: only scaffolds of size 1 can not be oriented
# Approach 2: looking if the observed homologous adjacencies all agree on an orientation
# c is a correction: 1 for first gene of an adjacency, -1 for second gene
def orient_scf(gene,instance,c):
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
               
    
# Importing DeClone adjacencies without orientation
DECLONE_INSTANCES_AUX={}
def DECLONE_import_aux():
    declone_file.sort(key=itemgetter(5))
    prev_instance=("","","")    
    for l in declone_file:
        (l10,l11,l12,l13,l14,l15)=split_adjacency(l)
        gene1=l14[1]
        gene2=l14[2]
        adj=((l10,l14[0]),(l11[0],l11[1],l14[1]),(l12[0],l12[1],l14[2]),l13,(l15[0],l15[1],l15[2]),("NA","NA","NA","NA","NA")) 
        tl15=tuple(l15)
        if tl15!=prev_instance:
            DECLONE_ID.append(tl15)
            DECLONE_INSTANCES_AUX[tl15]=[]
        DECLONE_INSTANCES_AUX[tl15].append(adj)
        prev_instance=tl15

# Orienting DeClone adjacencies
def DECLONE_orient_adjacencies():
    for instance in DECLONE_instances():
        DECLONE_INSTANCES[instance]=[]
        for adj in DECLONE_INSTANCES_AUX[instance]:
            gene1=DECLONE_adj_gene1_name(adj)
            gene2=DECLONE_adj_gene2_name(adj)
            if DECLONE_adj_score(adj)<1:
                or1=orient_scf(gene1,instance,1)
                or2=orient_scf(gene2,instance,-1)
            else:
                (or1,or2)=(0,0)
            DECLONE_INSTANCES[instance].append((adj[0],adj[1],adj[2],adj[3],adj[4],adj[5],(or1,or2)))


def DECLONE_import():
    DECLONE_import_aux()
    DECLONE_orient_adjacencies()

# --------------- Importing files ---------------------------------

def import_files():
    genes_file1=open(sys.argv[1],"r").readlines()
    for l in genes_file1:
        if l[0]!="#":
            l1=l.rstrip().split("\t")
            genes_file.append((l1[0],l1[1],l1[2],l1[3],l1[4],int(l1[5]),int(l1[6])))

    besst_file1=open(sys.argv[2],"r").readlines()
    for l in besst_file1:
        if l[0]!="#":
            besst_file.append(l.rstrip().split("\t"))

    declone_file1=open(sys.argv[3],"r").readlines()
    for l in declone_file1:
        if l[0]!="#":
            declone_file.append(l.rstrip().split("\t"))

    GENES_import()
    OG_import()
    GENOMES_import()
    BESST_import()
    DECLONE_import()

# --------------- MAIN ---------------------------------

import_files()

if sys.argv[4]=="-s":
    target_species=sys.argv[5]
    output_file=open(sys.argv[6],"w")
    if not target_species in GENOMES_list():
        print "Usage: gene_file BESST_file DeClone_file -s species output_file"
    else:
        output_file.write("#"+str(sys.argv)+"\n")
        output_file.write("#> species:ctg1:gene1:ctg2:gene2:orientation_ctg1:orientation_ctg2\n#species_id:species_name\tctg1:lg_ctg1:gene1:pos:orient.\tctg2:lg_ctg2:gene2:pos:orientat.\tscore\instance=tree1:tree2:root1-root2\tbesst=vscore:dscore:nb_reads:orient_ctg1:orient_ctg2:distance_ctgs:distance_genes\n")
        for instance in DECLONE_instances():
            for adj in DECLONE_adj_per_instance(instance):
                species1=DECLONE_adj_species_id(adj)
                species2=DECLONE_adj_species_name(adj)
                species=species1+":"+species2
                instance=DECLONE_adj_instance(adj)
                if DECLONE_adj_score(adj)<1.0 and species2==target_species:
                    prefix=adjacency_write_short(adj)
                    gene1=DECLONE_adj_gene1_name(adj)
                    gene2=DECLONE_adj_gene2_name(adj)
                    output_file.write("> "+prefix+"\n")
                    for adj2 in DECLONE_adj_per_instance(instance):
                        if DECLONE_adj_species_name(adj2)!="ANCESTRAL":
                            output_file.write("\t"+adjacency_write(adj2)+"\n")
