# Computing adjacencies evolution from DeClone results

import sys
from Bio import Phylo

# READING ARGUMENTS --------------------------------------------------------
gtrees_file=sys.argv[1]
adjacencies_file=sys.argv[2]
output_file=sys.argv[3]

# DATA STRUCTURES ----------------------------------------------------------

NODES={} # Indexed by (tree,node_ID), entries = (species,event,name)

def node_sp((t,c)):
    return NODES[(t,c)][0]
def node_ev((t,c)):
    return NODES[(t,c)][1]
def node_name((t,c)):
    return NODES[(t,c)][2]
def node_set((t,c),(s,e,n)):
    NODES[(t,c)]=(s,e,n)
def node_list():
    return NODES.keys()
def node_print1((t,c)):
    return t+":"+c

# Structure recording the speciation/extant descendants of each speciation node
DESCENDANTS={}  # Indexed by pairs (tree,node_ID), then each entry is indexed by a species_ID and is a list of pairs (tree,node_ID)
PARENTS={}      # Indexed by pairs (tree,node_ID), then each entry is a pair (tree,node_ID)

def desc_init((t,c),s1,s2):
    DESCENDANTS[(t,c)]={}
    DESCENDANTS[(t,c)][s1]=[]
    DESCENDANTS[(t,c)][s2]=[]
def desc_get_species((t,c)):
    return DESCENDANTS[(t,c)].keys()
def desc_get((t,c),s):
    return DESCENDANTS[(t,c)][s]
def desc_add((t,c),s,(t1,c1)):
    DESCENDANTS[(t,c)][s].append((t1,c1))
def parents_set((t,c),(t1,c1)):
    PARENTS[(t,c)]=(t1,c1)
def parents_get((t,c)):
    return PARENTS[(t,c)]

# Structure recording the adjacencies a gene belongs to
ADJACENCIES={} # Indexed by pairs (tree,node_ID) each entry is of the form (tree,node_ID,score,instance) 
               # 

def adj_tree(a):
    return a[0]
def adj_node(a):
    return a[1]
def adj_score(a):
    return a[2]
def adj_instance(a):
    return a[3]
def adj_init((t,c)):
    ADJACENCIES[(t,c)]=[]
def adj_add((t,c),(t1,c1,s,inst)):
    ADJACENCIES[(t,c)].append((t1,c1,s,inst))
def adj_get((t,c)):
    return ADJACENCIES[(t,c)]
def adj_list():
    return ADJACENCIES.keys()
def adj_print1(a):
    return a[0]+":"+a[1]+":"+str(a[2])+":"+a[3]


# BLOCK 1: Manipulation of Newick trees ------------------------------------

# Functions extracting information from Newick comments
def clade_ev(clade): # Event associated to a node
    c=clade.comment.split(":")
    for l in c:
        if l[0:3]=="Ev=":            
            return l[3:]
def clade_sp(clade): # Species associated to a node
    c=clade.comment.split(":")
    for l in c:
        if l[0:2]=="S=":            
            return l[2:]
def clade_nd(clade): # ID of a node
    c=clade.comment.split(":")
    for l in c:
        if l[0:3]=="ND=":            
            return l[3:]

# Functions to filter nodes in a tree traversal
def find_ext_spec(clade): # Keeping only extant or speciation nodes
    ev=clade_ev(clade)
    return ev=="Extant" or ev=="Spec"
def find_ext(clade):
    ev=clade_ev(clade)
    return ev=="Extant"
def find_spec(clade):
    ev=clade_ev(clade)
    return ev=="Spec"

# BLOCK 2: Looking for descendants --------------------------------------

# Functions to find the pre-speciations descendants of a pre-speciation node
def find_descendants_aux(tree, clade, origin):
    ev=clade_ev(clade)
    if ev=="Spec" or ev=="Extant":
        desc_add((tree,clade_nd(origin)),clade_sp(clade),(tree,clade_nd(clade)))
        parents_set((tree,clade_nd(clade)),(tree,clade_nd(origin)))
    elif ev=="GDup":
        find_descendants_aux(tree, clade[0],origin)
        find_descendants_aux(tree, clade[1],origin)
def find_descendants(tree,clade):
    # Assumption: clade is a speciation node
    find_descendants_aux(tree, clade[0], clade)
    find_descendants_aux(tree, clade[1], clade)

# BLOCK 3: populating the NODES/DESCENDANTS/PARENTS structures -----------
#          instanciating the ADJACENCIES structure

def read_gene_trees(gtrees_file):
    gtrees=Phylo.parse(gtrees_file,"newick")
    tree=0 # Tree ID
    for gtree in gtrees:
        tree_ID=str(tree)
        clades=gtree.find_clades(find_ext_spec)
        for c in clades:
            node_set((tree_ID,clade_nd(c)),(clade_sp(c),clade_ev(c),c.name))
            parents_set((tree_ID,clade_nd(c)),("",""))
            adj_init((tree_ID,clade_nd(c)))
            if clade_ev(c)=="Spec":
                desc_init((tree_ID,clade_nd(c)),clade_sp(c[0]),clade_sp(c[1]))
        tree+=1
    gtrees=Phylo.parse(gtrees_file,"newick")
    tree=0 # Tree ID
    for gtree in gtrees:
        tree_ID=str(tree)
        clades=gtree.find_clades(find_ext_spec)
        for c in clades:
            if clade_ev(c)=="Spec":
                find_descendants(tree_ID,c)
        tree+=1

# BLOCK 4: reading adjacencies -------------------------------------------

def extract_adj(adj):
    adj1=adj.rstrip().split("\t")
    (t1,c1)=(adj1[1].split(":")[0],adj1[1].split(":")[1])
    (t2,c2)=(adj1[2].split(":")[0],adj1[2].split(":")[1])
    if ((t1,c1)<(t2,c2)):
        return ((t1,c1),(t2,c2),float(adj1[3]),adj1[5])
    else:
        return ((t2,c2),(t1,c1),float(adj1[3]),adj1[5])

def read_adjacencies(adjacencies_file):
    adjacencies=open(adjacencies_file,"r").readlines()
    for adj in adjacencies:
        if adj[0]!="#":
            ((t1,c1),(t2,c2),score,instance)=extract_adj(adj)
            adj_add((t1,c1),(t2,c2,score,instance))
            
# BLOCK 5: checking adjacency evolution ----------------------------------

def check_adj_parent((t1,c1),(t2,c2),score,instance,output):
    (pt1,pc1)=parents_get((t1,c1))
    (pt2,pc2)=parents_get((t2,c2))
    s2=node_sp((t1,c1))
    if (pt1,pc1)>(pt2,pc2):
        (pt1,pc1),(pt2,pc2)=(pt2,pc2),(pt1,pc1)    
    if (pt1,pc1)!=("",""):
        s1=node_sp((pt1,pc1))
        l1_adj=adj_get((pt1,pc1))
        parent=node_print1((pt1,pc1))
        found=False
        for a in l1_adj:
            if (adj_tree(a),adj_node(a))==(pt2,pc2):
                found=True
                p_adj=adj_print1(a)                
        if not found:
            output.append((s1+":"+s2,"\tNA:NA:NA:NA:NA:NA--"+node_print1((t1,c1))+":"+adj_print1((t2,c2,score,instance))+"\n"))

def check_adj_desc((t1,c1),(t2,c2),score,instance,output):
    s1=node_sp((t1,c1))
    species=desc_get_species((t1,c1))
    desc1,desc2={},{}
    for s in species:
        output_str=""
        output_str1=s1+":"+s                
        desc1[s]=desc_get((t1,c1),s)
        desc2[s]=desc_get((t2,c2),s)
        found=False
        for (dt1,dc1) in desc1[s]:
            desc=node_print1((dt1,dc1))
            l1_adj=adj_get((dt1,dc1))
            for a in l1_adj:
                if (adj_tree(a),adj_node(a)) in desc2[s]:
                    found=True
                    output_str+="\t"+node_print1((t1,c1))+":"+adj_print1((t2,c2,score,instance))+"--"+desc+":"+adj_print1(a)
        if not found:
            output_str+=("\t"+node_print1((t1,c1))+":"+adj_print1((t2,c2,score,instance))+"--NA:NA:NA:NA:NA:NA\n")
        else:
            output_str+=("\n")
        output.append((output_str1,output_str))

def adj_evolution(adjacencies_file,output_file):
    adjacencies=open(adjacencies_file,"r").readlines()
    output_list=[]
    for adj in adjacencies:
        if adj[0]!="#":
            ((t1,c1),(t2,c2),score,instance)=extract_adj(adj)
            # Checking if the adjacency has a parent
            check_adj_parent((t1,c1),(t2,c2),score,instance,output_list)
            # Checking if the adjacency has descendants
            if node_ev((t1,c1))!="Extant":
                check_adj_desc((t1,c1),(t2,c2),score,instance,output_list)
    output_list.sort()
    output=open(output_file,"w")
    output.write("# Format: parent_species:descendant_species\tlist of pairs of parent--descendant adjacencies\n")
    output.write("#         adjacency format: tree1:node1:tree1:node2:score:instance\n")
    output.write("#         NA:NA:NA:NA:NA:NA means non parent or descendant\n")
    for o in output_list:
        output.write(o[0]+o[1])

# MAIN -------------------------------------------------------------------

read_gene_trees(gtrees_file)
read_adjacencies(adjacencies_file)
adj_evolution(adjacencies_file,output_file)
