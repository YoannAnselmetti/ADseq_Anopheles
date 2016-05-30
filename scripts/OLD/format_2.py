# Formatting the results to highlight the extant gap lengths

import sys

input_file=open(sys.argv[1],"r").readlines()
output_file=open(sys.argv[2],"w")

# --------------------- Reading the input file ----------------------
for l in input_file:
    if l[0]=="#" and l[1]=="\t":
        l1=l.rstrip().split("\t")
        for f in l1:
            if f!="#":
                f1=f.split(":")
                if f1[1]=="1.0":
                    output_file.write("\t1.0:"+f1[6])
                elif f1[6]!="NA":
                    output_file.write("\t"+f1[1]+":"+f1[6])
            else:
                output_file.write("#")
        output_file.write("\n")
    else:
        output_file.write(l)
