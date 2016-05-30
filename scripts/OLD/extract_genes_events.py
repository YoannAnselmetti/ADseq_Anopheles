import sys

input_file=open(sys.argv[1],"r").readlines()
output_file=open(sys.argv[2],"w")

genes_list=[]
for l in input_file:
    if l[0]=="G":
        tree=l.split()[2].split("|")[0]
    elif l[0]=="(":
        l1=l.replace("["," ").replace("]"," ").split()
        for c in l1:
            if c[0]=="&":
                c1=c.split(":")
                genes_list.append(c1[3].replace("S=","")+"\t"+tree+":"+c1[2].replace("ND=","")+"\t"+c1[1].replace("Ev=","")+"\n")

genes_list.sort()

for g in genes_list:
    output_file.write(g)
