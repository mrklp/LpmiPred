import numpy as np
import os
import re
import sys
from tr import tr
import random
from collections import Counter
from getopt import getopt

######################################## Parameters ########################################

usage = '''
##################################################################################################
#                                                                                                #
#  Usage: python de_barcode_adapter.py  -q file.fastq -m  miRNAs.txt -s species  -o clean.fastq  #
#                                                                                                #
##################################################################################################
'''

opts,args = getopt(sys.argv[1:],'q:m:s:o:h',['help'])
for opt_name, opt_value in opts:
    if opt_name in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt_name == "-s":
        species = opt_value
    elif opt_name == "-m":
        miR_file = opt_value
    elif opt_name == "-q":
        fastq_file = opt_value
    elif opt_name == "-o":
        out_file = opt_value
    else:
        print(usage)
        sys.exit()

######################################## def ########################################

def get_miR_probe(miR_file,species):
    file = open(miR_file,'r')
    mature_dict = {}
    for line in file:
        if re.search(r'^>',line):
            id = re.search(r'^>(\S+)',line).group(1)
            sequence = ""
            for line in file:
                if re.search(r'^>',line):
                    sequence = tr('U','T',sequence)
                    mature_dict[id] = sequence
                    id = re.search(r'^>(\S+)',line).group(1)
                    sequence = ""
                    continue
                sequence += line.strip().upper()
    
    sequence = tr('U','T',sequence)
    mature_dict[id] = sequence
    file.close()
    
    mature_dict_ = {}
    for miR in mature_dict:
        if species in miR:
            mature_dict_[miR] = mature_dict[miR]
    
    miRs = list(mature_dict_.keys())
    random.shuffle(miRs)
    miR_list = [mature_dict_[key] for key in miRs[:100]]
    return miR_list


def get_barcode_adapter(miR_list,fastq_file):
    reads_num = 0
    seqs_5p = []
    seqs_3p = []
    for miR_seq in miR_list:
        with os.popen("head -n 4000000 %s |rg -N -m 1000 %s" % (fastq_file,miR_seq), 'r') as grep_results:
            miR_reads = grep_results.readlines()
        reads_num += len(miR_reads)
        for miR_read in miR_reads:
            miR_read = miR_read.strip()
            result_5p = re.search('^([ATCGN]+)%s'%miR_seq, miR_read)
            if result_5p:
                seqs_5p.append(result_5p.group(1))
            result_3p = re.search('%s([ATCGN]+)$'%miR_seq, miR_read)
            if result_3p:
                seqs_3p.append(result_3p.group(1))
    
    if reads_num > 1000:
        # 检测barcode
        barcode = 'Null'
        if len(seqs_5p) > 0.5*reads_num:
            result = Counter(seqs_5p)
            count_dict = {}
            for key in result:
                if len(key) in count_dict:
                    count_dict[len(key)] += result[key]
                else:
                    count_dict[len(key)] = result[key]
            
            keys_list = list(count_dict.keys())
            keys_list_sort = sorted(keys_list,key = lambda i:i,reverse=True)
            for key in keys_list_sort:
                if count_dict[key] > 0.5*reads_num:
                    barcode = key
                    break
        
        # 检测adapter
        if len(seqs_3p) > 0.5*reads_num:
            lng_list = [1 for seq_3p in seqs_3p if len(seq_3p)>5]
            if sum(lng_list) < 0.5*reads_num:
                adapter = 'Null'
            else:
                adapter = os.popen("dnapi.py %s" % fastq_file).readlines()[0].strip()
        else:
            adapter = 'Null'
    else:
        barcode = 0
        adapter = 0
    
    return barcode,adapter

######################################## Main ########################################


def main():
    miR_list = get_miR_probe(miR_file,species)
    barcode,adapter = get_barcode_adapter(miR_list,fastq_file)
    srr = fastq_file.split('/')[-1].split('.')[0]
    print("%s\t%s\t%s"%(srr,barcode,adapter))
    
    if barcode and adapter:
        if barcode == 'Null' and adapter != 'Null':
            os.system("fastp -a %s -i %s -o  %s --length_required 18 --length_limit 25 -q 20 -u 20 -3 -W 4 -M 20 -y -Y 30 -w 16" % (adapter,fastq_file,out_file))
            print(fastq_file+'---1')
        elif barcode == 'Null' and adapter == 'Null':
            os.system("cp %s %s"%(fastq_file,out_file))
            print(fastq_file+'---2')
        elif barcode != 'Null' and adapter != 'Null':
            os.system("fastp -a %s -f %s -i %s -o %s --length_required 18 --length_limit 25 -q 20 -u 20 -3 -W 4 -M 20 -y -Y 30 -w 16" % (adapter,barcode,fastq_file,out_file))
            print(fastq_file+'---3')
        else:
            os.system("fastp -A -f %s -i %s -o %s --length_required 18 --length_limit 25 -q 20 -u 20 -3 -W 4 -M 20 -y -Y 30 -w 16" % (barcode,fastq_file,out_file))
            print(fastq_file+'---4')
        
        os.system("rm %s"%fastq_file)

if __name__=='__main__': 
    main()

###################################################################################################