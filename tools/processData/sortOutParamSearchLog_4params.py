#! /usr/bin/env python

# Purpose: summarise the parameter search log for a specified folder

import numpy as np
from os import path
from sys import argv, exit

try:
    rootDir = path.abspath(argv[1])
    num_of_nodes = int(argv[2])
except:
    print "Usage: ./sortOutParamSearchLog.py path num_of_jobs"
    exit(-1)

# node structure
node_list = range(1,num_of_nodes+1)
logV2_filePattern = "param_search_log_v2.dat"
logV3_filePattern = "param_search_log_v3.dat"

# pre-allocate space
logV2_combined = np.array([]).reshape(0, 22)
logV3_combined = np.array([]).reshape(0, 22)

# loop over all log files
for inode in node_list:
    #print "start for node%d"%inode
    logV2_file = path.join(rootDir, 'node%d'%inode, logV2_filePattern)
    logV2_data = np.loadtxt(logV2_file)
    logV2_combined = np.concatenate((logV2_combined, logV2_data), axis=0)

    logV3_file = path.join(rootDir, 'node%d'%inode, logV3_filePattern)
    logV3_data = np.loadtxt(logV3_file)
    logV3_combined = np.concatenate((logV3_combined, logV3_data), axis=0)
    #print 'node %d collected!'%inode
    if logV2_data.shape[0]!=8 or logV3_data.shape[0]!=8:
	print inode


# write to file
filename = path.join("param_search_log_v2_%s.dat"%rootDir.split('/')[-1])
np.savetxt(filename, logV2_combined, fmt="%18.6f", 
    delimiter="\t")

filename = path.join("param_search_log_v3_%s.dat"%rootDir.split('/')[-1])
np.savetxt(filename, logV3_combined, fmt="%18.6f", 
    delimiter="\t")

