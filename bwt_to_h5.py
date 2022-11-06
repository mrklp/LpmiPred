import h5py
import numpy as np
import os
import sys
from getopt import getopt

######################################## Parameters ########################################
usage = '''
#################################################################
#                                                               #
#  Usage: python bwt_to_h5.py -s signature.bwt -o signature.h5  #
#                                                               #
#################################################################
'''

thread_counts = 28
opts,args = getopt(sys.argv[1:],'s:o:h',['help'])
for opt_name, opt_value in opts:
    if opt_name in ("-h", "--help"):
        print(usage)
        sys.exit()
    elif opt_name == "-s":
        file_in = opt_value
    elif opt_name == "-o":
        file_out = opt_value
    else:
        print(usage)
        sys.exit()

######################################## Main ########################################

def main():
    hash_db = {}
    with open(file_in,'r') as bwt_file:
        for line in bwt_file:
            db = line.split('\t')[2]
            if db not in hash_db:
                hash_db[db] = line.strip()
            else:
                hash_db[db] += '#'+line.strip()
    
    with h5py.File(file_out,'w') as h5_file:
        for db in hash_db:
            h5_file[db] = hash_db[db]

if __name__=='__main__': 
    main()

###################################################################################################