# Labeling an adjacency evolution file by the gambiae chromosomes

import sys

label_file=open(sys.argv[1],"r").readlines()
adj_evol_file=open(sys.argv[2],"r").readlines()
output_file=open(sys.argv[3],"w")

status={}
for l in label_file:
    if l[0]!="#":
        l1=l.rstrip().split("\t")
        status[l1[0]]=l1[1]
        

for l in adj_evol_file:
    if l[0]=="#":
        output_file.write(l)
    else:
        l1=l.split("\t")
        adj1=l1[1].split("--")
        adj0=adj1[0]
        adj1=adj1[1]
        if adj0[0]=="N":
            adj=adj1
        else:
            adj=adj0
        tree0=adj.split(":")[0]
        tree1=adj.split(":")[2]
        output_file.write(status[tree0]+":"+status[tree1]+"\t"+l)
        
