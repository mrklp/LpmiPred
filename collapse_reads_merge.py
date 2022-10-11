import numpy as np
import os
import sys
from getopt import getopt



######################################## Parameters ########################################

usage = '''
######################################################################################################
#                                                                                                    #
#  Usage: python collapse_reads_merge.py -f ./collapse_reads_single/ -o collapse_reads_all.fasta     #
#                                                                                                    #
######################################################################################################
'''

opts,args = getopt(sys.argv[1:],'a:o:h',['help'])
for opt_name, opt_value in opts:
    if opt_name in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt_name == "-f":
        fasta_path = opt_value
    elif opt_name == "-o":
        out_file = opt_value
    else:
        print(usage)
        sys.exit()

######################################## main ########################################

def main():
    fasta_ids = []
    for dirpath,dirnames,filenames in os.walk(fasta_path):
        for filename in filenames:
            fasta_ids.append(os.path.join(dirpath,filename))
    
    fasta_dict = {}
    for fasta_id in fasta_ids:
        if os.path.getsize(fasta_id):
            with open(fasta_id,"r") as file_fasta:
                for line in file_fasta:
                    if re.search(r'^>',line):
                        count = re.search(r'^>\S+_x(\d+)$',line).group(1)
                        sequence = ""
                        for line in file_fasta:
                            if re.search(r'^>',line):
                                if sequence in fasta_dict.keys():
                                    fasta_dict[sequence] += int(count)
                                else:
                                    fasta_dict[sequence] = int(count)
                                count = re.search(r'^>\S+_x(\d+)$',line).group(1)
                                sequence = ""
                                continue
                            sequence += line.strip().upper()
                if sequence in fasta_dict.keys():
                    fasta_dict[sequence] += int(count)
                else:
                    fasta_dict[sequence] = int(count)
    
    dict_sorted = sorted(fasta_dict.items(), key = lambda kv:(kv[1], kv[0]),reverse=True)
    
    j = 0
    with open(out_file,"w")as write_file:
        for sequence,count in dict_sorted:
            j+=1
            write_file.write('>%s_x%s\n' % (j,count))
            write_file.write('%s\n' % sequence)

if __name__=='__main__': 
    main()

###################################################################################################