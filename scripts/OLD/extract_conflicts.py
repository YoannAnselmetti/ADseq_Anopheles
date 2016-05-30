# Separating sure adjacencies from dubious ones

import sys

adj_file=open(sys.argv[1],"r").readlines()
output1_file=open(sys.argv[2],"w")
output2_file=open(sys.argv[3],"w")

for l in adj_file:
    if l[0]=="#":
        c=1 # Do nothing
    else:
        l1=l.split("\t")
        conf=l1[7].split(":")
        or1=l1[3]
        or2=l1[4]
        if or1=="?" or or2=="?" or (int(conf[0])+int(conf[1])>0):
            output2_file.write(l)
        else: 
            output1_file.write(l)            
