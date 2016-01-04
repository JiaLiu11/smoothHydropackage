#! /usr/bin/env python

# Purpose: combine the sfactor log of each nodes

import numpy as np
from os import path
from sys import argv, exit

try:
    rootDir = path.abspath(argv[1])
    num_of_nodes = int(argv[2])
except:
    print "Usage: ./sortOutSfactorlog.py path num_of_jobs"
    exit(-1)

# node structure
node_list = range(1,num_of_nodes+1)
sfactorLog_filePattern = "sfactor_log.dat"

def unique_rows(a):
    """
    Pick out the unique rows
    """
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# pre-allocate space
log_combined = np.array([]).reshape(0, 7)
# loop over all log files
for inode in node_list:
    sfactorLog_file = path.join(rootDir, 'node%d'%inode, sfactorLog_filePattern)
    log_data = np.loadtxt(sfactorLog_file)
    log_data_unique = unique_rows(log_data)
    # lines_useful = np.abs(log_data_unique[:,-1]-totaldEdy_standard)<10
    # log_useful   = log_data_unique[lines_useful, :]
    log_useful = log_data_unique[log_data_unique[:,0]==1,:]
    if log_useful.shape[0]!=8:
	print inode
    log_combined = np.concatenate((log_combined, log_useful), axis=0)

# write to file
header_line = "Matching Time \t eta/s \t edec \t visbulknorm \t sfactor \t totaldNdy"
filename = path.join("sfactor_log_%s.dat"%rootDir.split('/')[-1])
np.savetxt(filename, log_combined[:,1::], fmt="%18.6f", 
    delimiter="\t")
