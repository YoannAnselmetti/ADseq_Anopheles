import sys
import os

print sys.argv

INSTANCES_LIST_FILENAME = sys.argv[1]
INSTANCES_DIRNAME =       sys.argv[2]
OUTPUT_DIRNAME =          sys.argv[3] # No / at the end
LOG_FILENAME =            sys.argv[3]+"_log"
MIN_SIZE =                int(sys.argv[4])
MAX_SIZE =                int(sys.argv[5])
KT =                      sys.argv[6]
GENOMES =                 sys.argv[7]

os.system('rm -f '+LOG_FILENAME)

input_file=open(INSTANCES_LIST_FILENAME,"r").readlines()
for line1 in input_file:
    header=line1.rstrip()
    adj_class=header.split()[2].replace("|","_")
    TREES_FILENAME=INSTANCES_DIRNAME+"/"+adj_class+"_trees"
    ADJ_FILENAME=  INSTANCES_DIRNAME+"/"+adj_class+"_adj"
    GENOMES_FILENAME = GENOMES
    input_trees=open(TREES_FILENAME,"r").readlines()
    g1 = ''
    g2 = ''
    for line2 in input_trees:
        if "NHX" in line2:
            if not g1:
                g1 = line2
                s1 = g1.count("NHX")
            elif not g2:
                g2 = line2               
                s2 = g2.count("NHX")                
                if min(s1,s2)<MIN_SIZE:
                    os.system('echo \"'+header+' '+str(max(s1,s2))+' '+str(min(s1,s2))+'  NOT PROCESSED (SMALL TREE)\" >> '+LOG_FILENAME)
                elif min(s1,s2)>MAX_SIZE:
                    os.system('echo \"'+header+' '+str(max(s1,s2))+' '+str(min(s1,s2))+'  NOT PROCESSED (LARGE TREE)\" >> '+LOG_FILENAME)        
                else:
                    OUTPUT_FILENAME=OUTPUT_DIRNAME+"/"+adj_class+"_scaff"
                    os.system('echo \"'+header+' '+str(s1)+' '+str(s2)+' PROCESSED\" >> '+LOG_FILENAME)
                    os.system('echo \"'+header+'\" > '+OUTPUT_FILENAME)
                    tmp1 = adj_class+"_tmp1.nhx"
                    fd1 = open(tmp1, 'w')
                    fd1.write(g1)
                    fd1.close()
                    tmp2 = adj_class+"_tmp2.nhx"
                    fd2 = open(tmp2, 'w')
                    fd2.write(g2)
                    fd2.close()
                    os.system('echo "SCAFFOLDING ENSEMBL PROBABILITIES kT ' + KT+ '" >>'+ OUTPUT_FILENAME)
                    os.system('../DeClone -t1 '+tmp1+' -t2 '+tmp2+' -a ' + ADJ_FILENAME + ' -gs '+GENOMES_FILENAME+' -sc 1.0 1.0 -i -kT '+KT+' >>'+OUTPUT_FILENAME)
                    os.system('rm -f '+adj_class+'_tmp?.nhx')
    sys.stdout.flush()
