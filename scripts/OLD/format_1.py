# Formatting the results in analysis/*_all file in the format of the assembly group.
# Filtering by minimum score and discarding the adjacencies that circularize a contig/scaffold
# Generates 3 files: a long format file, a short formal file and a conflicts file

import sys

input_file=open(sys.argv[1],"r").readlines()
genes_file1=open(sys.argv[2],"r").readlines()
output_file=open(sys.argv[3],"w")
output_file_short=open(sys.argv[3]+"_short","w")
conflict_file=open(sys.argv[3]+"_conflicts","w")
backbone_file=open(sys.argv[3]+"_backbone","w")
threshold=float(sys.argv[4])
nbpairs=float(sys.argv[5])

# --------------------- Reading the input file ----------------------
input_list=[]
prev_adj=""
SCF_ORIENT={}
SCF_LIST=[]
for l in input_file:
    if l[0]==">":
        if prev_adj!="":
            input_list.append(adjacency)
        prev_adj=l
        adjacency={}
        adjacency[0]=l.rstrip().split(" ")[1].split(":")
        adjacency[1]=[]
    elif l[0]!="#":
        l1=l.rstrip().split("\t")
        species=l1[1].split(":")[1]
        gene1=l1[2].split(":")
        gene2=l1[3].split(":")
        score=l1[4]
        besst=l1[6].split(":")
        adjacency[1].append((species,gene1,gene2,score,besst))
        SCF_ORIENT[(gene1[0],"+")]=0
        SCF_ORIENT[(gene1[0],"-")]=0
        SCF_ORIENT[(gene1[0],"?")]=0
        if not gene1[0] in SCF_LIST:
            SCF_LIST.append(gene1[0])
        SCF_ORIENT[(gene2[0],"+")]=0
        SCF_ORIENT[(gene2[0],"-")]=0
        SCF_ORIENT[(gene2[0],"?")]=0
        if not gene2[0] in SCF_LIST:
            SCF_LIST.append(gene2[0])
input_list.append(adjacency)

# --------------- Reading genes file ---------------------------------
# Genes file format
# Anopheles_albimanus  KB672286  GF0003303  AALB004105  -  461220  466589; TO DO: add exon information

genes_file=[]
for l in genes_file1:
    if l[0]!="#":
        l1=l.rstrip().split("\t")
        genes_file.append((l1[0],l1[1],l1[2],l1[3],l1[4],int(l1[5]),int(l1[6])))

GENES_ID=[]     # List of all gene names
GENES_INFO={}   # Info on each gene
GENES=[GENES_ID,GENES_INFO]

def GENES_list():
    return GENES[0]
def GENE_species(gene):
    return GENES[1][gene][0]
def GENE_scf(gene):
    return GENES[1][gene][1]
def GENE_sign(gene):
    return GENES[1][gene][3]
def GENE_start(gene):
    return int(GENES[1][gene][4])
def GENE_end(gene):
    return int(GENES[1][gene][5])
def GENE_og(gene):
    return GENES[1][gene][2]

def GENES_import():
    for (species,scf,og,gene,sign,start,end) in genes_file:
        GENES[0].append(gene)
        GENES[1][gene]=(species,scf,og,sign,start,end)

# --------------- Processing the input: reformatting -----------------
rev={"+":"-","-":"+","?":"?","NA":"NA"}

def adj_write1(s,c1,c2,g1,g2,or1,or2,score,besst,conf1,conf2):
    if besst[5]=="NA":
        d="?"
    else:
        d=besst[5]
    return s+"\t"+c1+"\t"+c2+"\t"+or1+"\t"+or2+"\t"+d+"\t"+score+"\t"+conf1+":"+conf2+"\t"+besst[0]+":"+besst[1]+":"+besst[2]+":"+besst[3]+":"+besst[4]

def adj_write2(s,g1,scf1,g2,scf2,score,besst):
    if scf1==scf2:
        if GENE_start(g1)<GENE_start(g2):
            dg=str(GENE_start(g2)-GENE_end(g1)+1)
        else:
            dg=str(GENE_start(g1)-GENE_end(g2)+1)
    else:
        dg=besst[6]
    if dg=="NA":
        dg="?"
    return s+"\t"+scf1+":"+g1+"\t"+scf2+":"+g2+":"+dg+"\t"+score+"\t"+besst[0]+":"+besst[1]+":"+besst[2]+":"+besst[3]+":"+besst[4]

def adj_write3(s,g1,scf1,g2,scf2,score,besst):
    if scf1==scf2:
        if GENE_start(g1)<GENE_start(g2):
            dg=str(GENE_start(g2)-GENE_end(g1)+1)
        else:
            dg=str(GENE_start(g1)-GENE_end(g2)+1)
    else:
        dg=besst[6]
    return s+":"+score+":"+scf1+":"+g1+":"+scf2+":"+g2+":"+dg+":"+besst[2]

GENES_import()
output_file.write("#species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tdistance\tscore\tconflict_ctg1:conflict_ctg2\tbesst_vscore:besst_dscore:nb_reads:orientation1:orientation2\n")
output_file_short.write("#Format\n")
output_file_short.write("#Funestus adjacencies: species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tdistance\tscore\tconflict_ctg1:conflict_ctg2\tbesst_vscore:besst_dscore:nb_reads:orientation1:orientation2\n")
output_file_short.write("#Supporting oberved adjacencies:    species:1.0:ctg1:gene1:ctg2:gene2:dist_genes:NA\n")
output_file_short.write("#Supporting inferreded adjacencies: species:score:ctg1:gene1:ctg2:gene2:dist_genes:number_of_supporting_paired-reads\n")

for adj in input_list:
    (species,ctg1,gene1,ctg2,gene2,or1,or2)=adj[0]
    for (s,g1,g2,score,besst) in adj[1]:
        besst0=besst[0]
        besst2=besst[2]
        if (g1[2],g2[2])==(gene1,gene2) and g1[0]!=g2[0] and (float(score)>=threshold or (besst0!="NA" and int(besst2)>=nbpairs)):
            SCF_ORIENT[(ctg1,or1)]+=1
            SCF_ORIENT[(ctg2,rev[or2])]+=1

#for scf in SCF_LIST:
#    scf1=SCF_ORIENT[(scf,"+")]
#    scf2=SCF_ORIENT[(scf,"-")]
#    scf3=SCF_ORIENT[(scf,"?")]
#    if scf1>1 or scf2>1 or scf1+scf2+scf3>2:
#        conflicts_file.write("CONFLICT\t"+scf+"\t"+str(scf1)+":"+str(scf2)+":"+str(scf3)+"\n")

nbconf=0;
nbadj=0
for adj in input_list:
    str2print1=""
    str2print2=""
    str2print3=""
    str2print2_short="#"
    str2print3_short=""
    (species,ctg1,gene1,ctg2,gene2,or1,or2)=adj[0]
    if ctg1>ctg2:
        ctg1,ctg2=ctg2,ctg1
        gene1,gene2=gene2,gene1
        or1,or2=rev[or2],rev[or1]
        reversed=True
    else:
        reversed=False

    (conf1,conf2)=("0","0")
    scf1=SCF_ORIENT[(ctg1,"+")]
    scf2=SCF_ORIENT[(ctg1,"-")]
    scf3=SCF_ORIENT[(ctg1,"?")]
    if scf1>1 or scf2>1 or scf1+scf2+scf3>2:
        conf1="1"
    scf1=SCF_ORIENT[(ctg2,"+")]
    scf2=SCF_ORIENT[(ctg2,"-")]
    scf3=SCF_ORIENT[(ctg2,"?")]
    if scf1>1 or scf2>1 or scf1+scf2+scf3>2:
        conf2="1"
        
    is_in_conflict=(conf1!="0" or conf2!="0")
    is_unoriented=(or1=="?" or or2=="?")
    is_kept=False
    for (s,g1,g2,score,besst) in adj[1]:
        besst1=besst
        if reversed:
            g1,g2=g2,g1
            besst1=(besst[0],besst[1],besst[2],rev[besst[4]],rev[besst[3]],besst[5],besst[6])
        besst0=besst1[0]
        besst2=besst1[2]
        if (g1[2],g2[2])==(gene1,gene2) and g1[0]!=g2[0] and (float(score)>=threshold or (besst0!="NA" and int(besst2)>=nbpairs)):
            str2print1=adj_write1(species,ctg1,ctg2,gene1,gene2,or1,or2,score,besst1,conf1,conf2)+"\n"
            is_kept=True
        if (g1[2],g2[2])!=(gene1,gene2) and float(score)==1:
            str2print2+="#"+adj_write2(s,g1[2],g1[0],g2[2],g2[0],score,besst1)+"\n"
            str2print2_short+="\t"+adj_write3(s,g1[2],g1[0],g2[2],g2[0],score,besst1)
        elif (g1[2],g2[2])!=(gene1,gene2) and (float(score)>=threshold or (besst0!="NA" and int(besst2)>=nbpairs)):
            str2print3+="#"+adj_write2(s,g1[2],g1[0],g2[2],g2[0],score,besst1)+"\n"
            str2print3_short+="\t"+adj_write3(s,g1[2],g1[0],g2[2],g2[0],score,besst1)
    if is_kept:
        output_file.write(str2print1+str2print2+str2print3)
        output_file_short.write(str2print1+str2print2_short+str2print3_short+"\n")

    if is_kept and not (is_in_conflict or is_unoriented):
        backbone_file.write(str2print1)
        nbadj+=1
    elif is_kept:
        conflict_file.write(str2print1)
        nbconf+=1

print str(threshold)+":\t"+str(nbadj+nbconf)+" adjacencies\twith "+str(nbconf)+" conflicting or unoriented"

#['Anopheles_funestus', 'KB668996', 'AFUN009019', 'KB668756', 'AFUN001925', '-', '+']



                                
