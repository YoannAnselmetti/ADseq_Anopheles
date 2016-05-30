# Filtering for pre-speciation adjacencies and reformating the DeClone output

import sys

declone_file=open(sys.argv[1],"r").readlines()
pre_speciation_file=open(sys.argv[2],"r").readlines()
output_file=open(sys.argv[3],"w")

genes_events={}
for l in pre_speciation_file:
    l1=l.rstrip().split()
    genes_events[(l1[0],l1[1])]=l1[2]

adjacencies=[]
for l in declone_file:
  l1=l.rstrip().split()
  if l1[0]=="Adjacency":
      class_id=l1[2]
      class_id1=class_id.split("|")
      tree1=class_id1[0]
      tree2=class_id1[1]
  elif l1[0]!="SCAFFOLDING" and l1[0]!="#T1:" and l1[0]!="#T2:":
      l2=l1[0].split(":")
      event1=genes_events[(l2[0],tree1+":"+l2[1])]
      event2=genes_events[(l2[0],tree2+":"+l2[2])]
      if (event1=="Extant" and event2=="Extant") or (event1=="Spec" and event2=="Spec"):
          if len(l1)==1:
              suffix="ANCESTRAL:NA:NA\t"+class_id+"\n"
          else:
              suffix=l1[1]+"\t"+class_id+"\n"
          adjacencies.append(l2[0]+"\t"+tree1+":"+l2[1]+"\t"+tree2+":"+l2[2]+"\t"+l2[3]+"\t"+suffix)
          #adjacencies.append(l2[0]+"\t"+tree2+":"+l2[2]+"\t"+tree1+":"+l2[1]+"\t"+l2[3]+"\t"+suffix)
      
adjacencies.sort()

output_file.write("#species ID\ttree1:node1\ttree2:node2\tDeClone score\tspecies:gene1:gene2\tinstance=tree1|tree2|root1-root2\n")
for a in adjacencies:
    output_file.write(a)
