#! /usr/bin/env python

import sys, shutil
from numpy import *
from os import path, makedirs, remove
from subprocess import call

def checkAvgMultiplicityandTransverseProfile():
   runRecord = open(path.abspath('./runRecord.dat'), 'a')
   errRecord = open(path.abspath('./errRecord.dat'), 'a')
   cmd = './superMC.e'
   nev = 1000
   sdsum = 0.0
   sdsum2 = 0.0
   sdavg = zeros([261, 261])
   for iev in range(nev):
     print("event : %d" % iev)
     if path.isfile(path.abspath('./data/sd_event_1_block.dat')) :
        remove(path.abspath('./data/sd_event_1_block.dat'))
     call(cmd, shell=True, stdout = runRecord, stderr = errRecord)
     sd = loadtxt(path.abspath('./data/sd_event_1_block.dat'))
     sdsum = sdsum + sum(sd)
     sdsum2 = sdsum2 + sum(sd)*sum(sd)
     sdavg = sdavg + sd
   meansd = sdsum/nev
   stdsd = sqrt(sdsum2/nev - meansd**2)
   print(meansd, stdsd)
   savetxt(path.abspath('./data/sd_event_avg_check.dat'), sdavg/nev, fmt='%.18e', delimiter='   ')


if __name__ == "__main__":
   checkAvgMultiplicityandTransverseProfile()
