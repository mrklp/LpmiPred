#!/usr/bin/python

import re
import sys
import os
from tr import tr
from getopt import getopt


######################################## Parameters ########################################
usage = '''
############################################################################################
#                                                                                          #
#  Usage: python excise_precursors.py -a 20  -b reads_mappings.bwt -g genome.fa -c chr1    #
#                                                                                          #
############################################################################################
'''

freq_min = 20
match,umatch,unmatch,space = [8,1,-3,-10]
count_excisions=0

opts,args = getopt(sys.argv[1:],'a:b:g:h',['help'])
for opt_name, opt_value in opts:
    if opt_name in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt_name == "-a":
        freq_min = int(opt_value)
    elif opt_name == "-b":
        file_bwt = opt_value
    elif opt_name == "-g":
        file_fasta = opt_value
    elif opt_name == "-c":
        chr = opt_value
    else:
        print(usage)
        sys.exit()


######################################## Main ########################################

hash_pos = {}
hash_genome = {}
parse_file_bwt(file_bwt)
parse_genome(file_fasta)
excise_precursors(species,chr)


######################################## Functions ########################################

def SmithWaterman(seq_1,seq_2,match,umatch,unmatch,space):
    seq_2 = seq_2[seq_2.index(seq_1)+len(seq_1):] if seq_2.index(seq_1)<=10 else seq_2[:seq_2.index(seq_1)]
    seq_1 = ['AG' if i== 'T' else 'CT' if i == 'G' else 'G' if i == 'C' else 'T' for i in list(seq_1[::-1])]
    row_num = len(seq_1)+2
    col_num = len(seq_2)+2
    
    matrix_score = [[0 for j in range(col_num)] for i in range(row_num)]
    matrix_arrow = [[0 for j in range(col_num)] for i in range(row_num)]
    
    for i in [0,1]:
        for j in range(2,col_num):
                matrix_score[i][j] = seq_2[j-2] if i==0 else 0
                matrix_arrow[i][j] = seq_2[j-2] if i==0 else 2
    
    for j in [0,1]:
        for i in range(2,row_num):
                matrix_score[i][j] = seq_1[i-2] if j==0 else 0
                matrix_arrow[i][j] = seq_1[i-2] if j==0 else 2
    
    for i in range(2,row_num):
        for j in range(2,col_num):
            road_1_score = matrix_score[i][j-1] + space
            road_2_score = matrix_score[i-1][j] + space
            if matrix_score[0][j] in matrix_score[i][0]:
                if matrix_score[i][0].index(matrix_score[0][j])==0:
                    road_3_score = matrix_score[i-1][j-1] + match
                if matrix_score[i][0].index(matrix_score[0][j])==1:
                    road_3_score = matrix_score[i-1][j-1] + umatch
            else:
                road_3_score = matrix_score[i-1][j-1] + unmatch
            if road_1_score > road_2_score:
                matrix_score[i][j] = road_1_score if road_1_score > road_3_score else road_3_score
            else:
                matrix_score[i][j] = road_2_score if road_2_score > road_3_score else road_3_score
            matrix_arrow[i][j] = 0 if matrix_score[i][j] == road_3_score else (1 if matrix_score[i][j] == road_2_score else -1)
            if matrix_score[i][j] < 0:
                matrix_score[i][j] = 0
                matrix_arrow[i][j] = 2
    
    #get a maximum alignment score
    row,col = [1,1]
    for i in range(2,row_num):
        for j in range(2,col_num):
            if matrix_score[i][j] > matrix_score[row][col]:
                row,col = [i,j] 
    
    max_align_score = matrix_score[row][col]
    
    row_,col_ = [1,1]
    for i in range(2,row_num):
        for j in range(2,col_num):
            if matrix_score[row][col] > matrix_score[i][j] > matrix_score[row_][col_] and (j<col-20 or j>col+20):
                row_,col_ = [i,j]
    
    sub_align_score = matrix_score[row_][col_]
    
    #get two aligned sequences
    align_seq_1 = ''
    align_num = 0
    align_seq_2 = ''
    i,j = [row,col]
    while (matrix_score[i][j]):
        if matrix_arrow[i][j] == -1:
            align_seq_1 += '-'
            align_seq_2 += matrix_arrow[0][j]
            j-=1
        elif matrix_arrow[i][j] == 1:
            align_seq_1 += matrix_arrow[i][0]
            align_num += 1
            align_seq_2 += '-'
            i-=1
        else:
            align_num += 1
            align_seq_1 += matrix_arrow[i][0]
            align_seq_2 += matrix_arrow[0][j]
            i-=1
            j-=1
    
    align_seq_1_ = ''
    align_num_ = 0
    align_seq_2_ = ''
    i,j = [row_,col_]
    while (matrix_score[i][j]):
        if matrix_arrow[i][j] == -1:
            align_seq_1_ += '-'
            align_seq_2_ += matrix_arrow[0][j]
            j-=1
        elif matrix_arrow[i][j] == 1:
            align_seq_1_ += matrix_arrow[i][0]
            align_num_ += 1
            align_seq_2_ += '-'
            i-=1
        else:
            align_num_ += 1
            align_seq_1_ += matrix_arrow[i][0]
            align_seq_2_ += matrix_arrow[0][j]
            i-=1
            j-=1

    align_seq_1 = align_seq_1[::-1]
    align_seq_2 = align_seq_2[::-1]
    align_seq_1_ = align_seq_1_[::-1]
    align_seq_2_ = align_seq_2_[::-1]
    
    return align_seq_1,align_seq_2,align_seq_1_,align_seq_2_,align_num,align_num_

def parse_file_bwt(file_bwt):
    file = open(file_bwt,'r')
    for line in file:
        line = line.strip().split('\t')
        query = line[0]
        db = line[2]
        strand = line[1]
        db_beg = str(int(line[3])+1)
        db_end = str(int(line[3])+len(line[4]))
        freq = int(re.search(r'_x(\S+)$',query).group(1))
        if db not in hash_pos:
            hash_pos.update({db:{strand:{db_beg:{db_end:freq}}}})
        elif strand not in hash_pos[db]:
            hash_pos[db].update({strand:{db_beg:{db_end:freq}}})
        elif db_beg not in hash_pos[db][strand]:
            hash_pos[db][strand].update({db_beg:{db_end:freq}})
        elif db_end not in hash_pos[db][strand][db_beg]:
            hash_pos[db][strand][db_beg].update({db_end:freq})
        else:
            hash_pos[db][strand][db_beg][db_end]+=freq
    file.close()
    return

def find_freq_max_downstream(db,strand,db_beg,db_end):
    db_beg=int(db_beg)
    db_end=int(db_end)
    freq_max=0
    pos_beg=db_beg+1
    while pos_beg<=db_end+70:
        if str(pos_beg) in hash_pos[db][strand].keys():
            for pos_end in [str(i) for i in sorted([int(j) for j in hash_pos[db][strand][str(pos_beg)].keys()])]:
                freq=hash_pos[db][strand][str(pos_beg)][str(pos_end)]
                if freq>freq_max:
                    freq_max=freq
        pos_beg+=1
    return freq_max

def outprecursor_(db,strand,excise_beg,excise_end,seq_1,seq_2,align_seq_2,direc,file_1,file_2):
    align_seq_2=align_seq_2.replace('-','')
    align_beg=seq_2.index(align_seq_2)+1
    align_end=align_beg+len(align_seq_2)-1
    global count_excisions
    if direc=='down':
        if(len(seq_2)-align_end>=10):
            excise_seq=seq_2[:align_end+10]
            if strand=="-":
                excise_seq = tr('acgtuACGTU','TGCAATGCAA',excise_seq)[::-1]
            new_end=excise_beg+align_end+10-1;
            count_excisions+=1
            file_1.write('>%s_%d\n%s\n'%(db,count_excisions,excise_seq))
            file_2.write('>%s_%d\t%s\t%d\t%d\n'%(db,count_excisions,strand,excise_beg,new_end))
        else:
            excise_seq=seq_2
            if strand=="-":
                excise_seq = tr('acgtuACGTU','TGCAATGCAA',excise_seq)[::-1]
            count_excisions+=1
            file_1.write('>%s_%d\n%s\n'%(db,count_excisions,excise_seq))
            file_2.write('>%s_%d\t%s\t%d\t%d\n'%(db,count_excisions,strand,excise_beg,excise_end))
    if direc=='up':
        if(align_beg<=10):
            excise_seq=seq_2
            if strand=="-":
                excise_seq = tr('acgtuACGTU','TGCAATGCAA',excise_seq)[::-1]
            count_excisions+=1
            file_1.write('>%s_%d\n%s\n'%(db,count_excisions,excise_seq))
            file_2.write('>%s_%d\t%s\t%d\t%d\n'%(db,count_excisions,strand,excise_beg,excise_end))
        else:
            excise_seq=seq_2[align_beg-10-1:]
            if strand=="-":
                excise_seq = tr('acgtuACGTU','TGCAATGCAA',excise_seq)[::-1]
            new_beg=excise_beg+align_beg-10-1
            count_excisions+=1
            file_1.write('>%s_%d\n%s\n'%(db,count_excisions,excise_seq))
            file_2.write('>%s_%d\t%s\t%d\t%d\n'%(db,count_excisions,strand,new_beg,excise_end))
    return

def outprecursor(db,strand,excise_beg,excise_end,seq_1, seq_2, direc,file_1,file_2):
    align_seq_1,align_seq_2,align_seq_1_,align_seq_2_,align_num,align_num_ = SmithWaterman(seq_1,seq_2,match,umatch,unmatch,space)
    align_ratio = align_num/len(seq_1)
    align_ratio_ = align_num_/len(seq_1)
    diff_lng=len(seq_1)-len(align_seq_2.replace('-',''))
    diff_lng_=len(seq_1)-len(align_seq_2_.replace('-',''))
    if align_ratio > 0.6 and -6<=diff_lng<=6:
        outprecursor_(db,strand,excise_beg,excise_end,seq_1,seq_2,align_seq_2,direc,file_1,file_2)
    if align_ratio_ > 0.6 and -6<=diff_lng_<=6:
        outprecursor_(db,strand,excise_beg,excise_end,seq_1,seq_2,align_seq_2_,direc,file_1,file_2)
    return

def print_positions(db,strand,db_seq,db_beg,db_end,file_1,file_2):
    excise_beg=max(1,int(db_beg)-10)
    excise_end=min(int(db_end)+1000,len(db_seq))
    seq_1=db_seq[int(db_beg)-1:int(db_end)]
    seq_2=db_seq[excise_beg-1:excise_end]
    if not (re.search('[^ACGTU]',seq_1) or re.search('[^ACGTU]',seq_2)):
        outprecursor(db,strand,excise_beg,excise_end,seq_1,seq_2, 'down',file_1,file_2)
    
    excise_beg=max(1,int(db_beg)-1000)
    excise_end=min(int(db_end)+10,len(db_seq))
    seq_1=db_seq[int(db_beg)-1:int(db_end)]
    seq_2=db_seq[excise_beg-1:excise_end]
    if not (re.search('[^ACGTU]',seq_1) or re.search('[^ACGTU]',seq_2)):
        outprecursor(db,strand,excise_beg,excise_end,seq_1,seq_2, 'up',file_1,file_2)
    return

def excise(db,db_seq,file_1,file_2):
    if db in hash_pos:
        for strand in sorted(hash_pos[db]):
            db_limit=0
            for db_beg in [str(i) for i in sorted([int(j) for j in hash_pos[db][strand].keys()])]:
                for db_end in [str(i) for i in sorted([int(j) for j in hash_pos[db][strand][db_beg].keys()])]:
                    freq=hash_pos[db][strand][db_beg][db_end]
                    freq_max_ds=find_freq_max_downstream(db,strand,db_beg,db_end)
                    if freq<freq_min or freq<freq_max_ds or int(db_beg)<db_limit:
                        continue
                    print_positions(db,strand,db_seq,db_beg,db_end,file_1,file_2)
                    db_limit=int(db_end)+70
    return

def parse_genome(file_fasta):
    file = open(file_fasta,'r')
    for line in file:
        if re.search(r'^>',line):
            id = re.search(r'^>(\S+)',line).group(1)
            sequence = ""
            for line in file:
                if re.search(r'^>',line):
                    hash_genome[id] = sequence
                    id = re.search(r'^>(\S+)',line).group(1)
                    sequence = ""
                    continue
                sequence += line.strip().upper()
    hash_genome[id] = sequence
    file.close()
    return

def excise_precursors(species,db):
    file_1 = open('precursors_%s.fa'%db,'w')
    file_2 = open('precursors_%s.coords'%db,'w')
    excise(db,hash_genome[db],file_1,file_2)
    file_1.close()
    file_2.close()
    return

###################################################################################################