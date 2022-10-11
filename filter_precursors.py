#!/usr/bin/python

import os
import re
import sys
import h5py
from tr import tr
from math import ceil
from getopt import getopt

######################################## Parameters ########################################
usage = '''
###########################################################################################################################################
#                                                                                                                                         #
#  Usage: python 13_filter_precursors.py  -n 20 -f precursors.fa -c precursors.coords  -s signature.h5 -p precursors.str  -o result.txt   #
#                                                                                                                                         #
###########################################################################################################################################
'''

high_conf = 20
opts,args = getopt(sys.argv[1:],'n:s:t:i:h',['help'])
for opt_name, opt_value in opts:
    if opt_name in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt_name == "-n":
        high_conf = int(opt_value)
    elif opt_name == "-f":
        file_precursors = opt_value
    elif opt_name == "-c":
        file_coords = opt_value
    elif opt_name == "-s":
        file_signature = opt_value
    elif opt_name == "-p":
        file_structure = opt_value
    elif opt_name == "-o":
        file_result = opt_value
    else:
        print(usage)
        sys.exit()

######################################## Functions ########################################

def parse_file_coords(file_coords):
    hash_coords = {}
    file_in = open(file_coords,'r')
    for line in file_in:
        query,strand,beg,end = line.strip().split('\t')
        chr = re.search(r'^>(\S+)_',line).group(1)
        hash_coords[query[1:]]=[chr,strand,int(beg),int(end)]
    file_in.close()
    return  hash_coords

def parse_file_struct(file_structure):
    hash_seq = {}
    hash_struct = {}
    hash_mfe = {}
    file_in = open(file_structure,'r')
    for line in file_in:
        line = line.strip()
        if re.search(r'^>',line):
            id = re.search(r'^>(\S+)',line).group(1)
            seq = ''
            struct = ''
            mfe = ''
            for line in file_in:
                line = line.strip()
                if re.search(r'^>',line):
                    hash_seq[id] = seq
                    hash_struct[id] = struct
                    hash_mfe[id] = mfe
                    id = re.search(r'^>(\S+)',line).group(1)
                    seq = ''
                    struct = ''
                    mfe = ''
                    continue
                if re.fullmatch('[ACGTU]+',line):
                    seq += tr('uU','tT',line)
                if re.search('([.()]+)',line):
                    struct += re.search('([.()]+)',line).group(1)
                if re.search('\(\s*(-*\d+.\d+)\)',line):
                    mfe = float(re.search('\(\s*(-*\d+.\d+)\)',line).group(1))
    hash_seq[id] = seq
    hash_struct[id] = struct
    hash_mfe[id] = mfe
    file_in.close()
    return  hash_seq,hash_struct,hash_mfe

def pass_filtering_structure(db,hash_bp,hash_comp):
    
    # potential mature and star must be identifiable
    if "p5_struct" in hash_comp and "p3_struct" in hash_comp:
        if not (re.fullmatch('[(.]+',hash_comp["p5_struct"]) and re.fullmatch('[.)]+',hash_comp["p3_struct"])):
            return
    else:
        return
    
    p5_beg = hash_comp["p5_beg"]
    p5_end = hash_comp["p5_end"]
    p3_beg = hash_comp["p3_beg"]
    p3_end = hash_comp["p3_end"]
    
    # minimum 60% base pairings in duplex
    pair_counts = 0
    for pos in range(p5_beg,p5_end+1):
        if str(pos) in hash_bp:
            if hash_bp[str(pos)] in range(p3_beg,p3_end+1):
                pair_counts += 1
    if not (pair_counts/len(hash_comp["p5_seq"])>=0.6 and pair_counts/len(hash_comp["p3_seq"])>=0.6):
        return
    
    # must with 0-4 nt overhang at 3' end
    for pos in range(p5_beg,p5_end+1):
        if str(pos) in hash_bp:
            p5_end_bp= pos
    for pos in range(p3_beg,p3_end+1):
        if str(pos) in hash_bp:
            p3_end_bp= pos
    if not (p5_end-p5_end_bp<=4 and p3_end-p3_end_bp<=4):
        return
    
    # must have a folding free energy of <-0.2kal/mol/nt
    if not (hash_mfe[db]/len(hash_struct[db])<-0.2):
        return
    
    return 1

def pass_filtering_signature(db,hash_query,hash_comp):
    p5_beg_freq = 0
    for query in hash_query:
        if hash_query[query]["db_beg"] == hash_comp["p5_beg"]:
            p5_beg_freq+=hash_query[query]["freq"]
    p3_beg_freq = 0
    for query in hash_query:
        if hash_query[query]["db_beg"] == hash_comp["p3_beg"]:
            p3_beg_freq+=hash_query[query]["freq"]
    
    if p5_beg_freq<=high_conf or p3_beg_freq<=high_conf:
        return 0,0
    db_lng = len(hash_seq[db])
    p5_all_freq = 0
    for query in hash_query:
        if hash_query[query]["db_beg"]>=max(1,hash_comp["p5_beg"]-20) and hash_query[query]["db_end"]<=min(hash_comp["p5_end"]+20,db_lng):
            p5_all_freq += hash_query[query]["freq"]
    
    p3_all_freq = 0
    for query in hash_query:
        if hash_query[query]["db_beg"]>=max(1,hash_comp["p3_beg"]-20) and hash_query[query]["db_end"]<=min(hash_comp["p3_end"]+20,db_lng):
            p3_all_freq += hash_query[query]["freq"]
    
    if not (p5_beg_freq/p5_all_freq >= 0.5 and p3_beg_freq/p3_all_freq >= 0.5):
        return 0,0
    
    return p5_beg_freq,p3_beg_freq

def fill_structure(db):
    hash_bp={}
    struct=hash_struct[db]
    lng=len(struct)
    bps=[]
    pos = 0
    for struct_pos in struct:
        pos+=1
        if struct_pos=='(':
            bps.append(pos)
        if struct_pos==')':
            pos_prev=bps.pop()
            hash_bp[str(pos_prev)]=pos
            hash_bp[str(pos)]=pos_prev
    return  hash_bp

def fill_pri(db):
    hash_comp={}
    seq=hash_seq[db]
    struct=hash_struct[db]
    mfe=hash_mfe[db]
    length=len(seq)
    
    hash_comp["pri_id"]=db
    hash_comp["pri_seq"]=seq
    hash_comp["pri_struct"]=struct
    hash_comp["pri_mfe"]=mfe
    hash_comp["pri_beg"]=1
    hash_comp["pri_end"]=length
    return hash_comp

def fill_miRNA(db,hash_query):
    hash_comp = {}
    p5_query,p3_query = find_miRNA_query(db,hash_query)
    if not (p5_query and p3_query):
        return hash_comp
    p5_beg = hash_query[p5_query]["db_beg"]
    p5_end = hash_query[p5_query]["db_end"]
    p5_strand = hash_query[p5_query]["strand"]
    p5_seq = hash_seq[db][p5_beg-1:p5_end]
    p5_struct = hash_struct[db][p5_beg-1:p5_end]
    
    p3_beg = hash_query[p3_query]["db_beg"]
    p3_end = hash_query[p3_query]["db_end"]
    p3_strand = hash_query[p3_query]["strand"]
    p3_seq = hash_seq[db][p3_beg-1:p3_end]
    p3_struct = hash_struct[db][p3_beg-1:p3_end]
    
    hash_comp["p5_query"] = p5_query
    hash_comp["p5_beg"] = p5_beg
    hash_comp["p5_end"] = p5_end
    hash_comp["p5_strand"] = p5_strand
    hash_comp["p5_seq"] = p5_seq
    hash_comp["p5_struct"] = p5_struct
    
    hash_comp["p3_query"] = p3_query
    hash_comp["p3_beg"] = p3_beg
    hash_comp["p3_end"] = p3_end
    hash_comp["p3_strand"] = p3_strand
    hash_comp["p3_seq"] = p3_seq
    hash_comp["p3_struct"] = p3_struct
    
    return  hash_comp

def find_miRNA_query(db,hash_query):
    db_lng = len(hash_seq[db])
    query_freq = {}
    for query in hash_query:
        freq = hash_query[query]["freq"]
        query_freq[query] = freq
    query_freq_sorted = sorted(query_freq.items(),key=lambda i:i[1],reverse=True)
    queries = [i[0] for i in query_freq_sorted]
    for query in queries:
        if hash_query[query]["db_end"] <=35:
            p5_query = query
            break
        else:
            p5_query = 0
    for query in queries:
        if db_lng-hash_query[query]["db_beg"]<=35:
            p3_query = query
            break
        else:
            p3_query = 0
    
    return  p5_query,p3_query

######################################## Main ########################################


def main():
    hash_coords = parse_file_coords(file_coords)
    hash_seq,hash_struct,hash_mfe = parse_file_struct(file_structure)
    
    precursors_id = os.popen("grep '>' %s"%file_precursors).readlines()
    precursors_id = [id.strip()[1:] for id in precursors_id]
    
    file_out = open(file_result,'w')
    h5_file = h5py.File(file_signature,'r')
    
    for db in precursors_id:
        
        hash_query={}
        db_query = []
        for key in h5_file.keys():
            if db in h5_file[key]:
                db_query += str(h5_file[key][db][...]).split("'")[1].split('#')
        
        if db_query == []:
            print('##%s##'%db)
            continue
        
        for line in db_query:
            query,strand,db,db_beg,db_seq = line.split('\\t')[:5]
            db_beg = int(db_beg)+1
            db_end = db_beg+len(db_seq)-1
            freq = int(re.search('_x(\d+)',query).group(1))
            hash_query[query]={"db_beg":int(db_beg),
                                "db_end":int(db_end),
                                "strand":strand,
                                "freq":int(freq)}
        hash_bp = fill_structure(db)
        hash_comp = fill_miRNA(db,hash_query)
        if pass_filtering_structure(db,hash_bp,hash_comp):
            p5_beg_freq,p3_beg_freq = pass_filtering_signature(db,hash_query,hash_comp)
            if p5_beg_freq!=0:
                chr,strand,beg,end = hash_coords[db]
                file_out.write('%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\n'%(db,hash_comp["p5_seq"],p5_beg_freq,hash_comp["p3_seq"],p3_beg_freq,hash_seq[db],chr,strand,beg,end))
        
        print('##%s##'%db)
    
    h5_file.close()
    file_out.close()

if __name__=='__main__': 
    main()

###################################################################################################



