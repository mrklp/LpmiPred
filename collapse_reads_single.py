import numpy as np
import os
import sys
from tr import tr
from collections import Counter
from getopt import getopt



######################################## Parameters ########################################

usage = '''
#####################################################################################
#                                                                                   #
#  Usage: python collapse_reads_single.py -a clean.fasta -o collapse_reads.fasta    #
#                                                                                   #
#####################################################################################
'''

opts,args = getopt(sys.argv[1:],'a:o:h',['help'])
for opt_name, opt_value in opts:
    if opt_name in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt_name == "-a":
        fasta_file = opt_value
    elif opt_name == "-o":
        out_file = opt_value
    else:
        print(usage)
        sys.exit()

######################################## main ########################################

def main():
    srx = fasta_file.split('/')[-1].split('.')[0]
    with os.popen("egrep '^[ATCGUNatcgun]+$' %s" % fasta_file) as grep_reads:
        reads = grep_reads.readlines()
    
    reads = [read.strip() for read in reads]
    reads_count = Counter(reads)
    reads_sort = sorted(reads_count.items(), key = lambda kv:(kv[1], kv[0]),reverse=True)
    with open(out_file,"w")as write_file:
        i = 0
        for read,count in reads_sort:
            if 18<=len(read)<=25 and count >1:
                i+=1
                read = tr('acgtunU','ACGTTNT',read)
                write_file.write('>%s_%s_x%s\n' % (srx,i,count))
                write_file.write('%s\n' % read)
        
    os.system("rm %s"%fasta_file)

if __name__=='__main__': 
    main()

###################################################################################################