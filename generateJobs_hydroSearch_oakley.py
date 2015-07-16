#! /usr/bin/env python
"""
    This script duplicates the EBE-Node folder and generate a collection of pbs
    files to be batch-submitted. For efficiency all codes inside EBE-Node should
    be compiled.
"""

from sys import argv, exit
from os import makedirs, path, unlink
from shutil import copytree, copy, rmtree

from check_prerequisites import check_environment, check_executables, greetings

# check argv
estimatedRunTime = 167 # oakley allows 168 Hrs for serial jobs
try:
    # set parameters
    numberOfJobs = int(argv[1])
    # set optional parameters
    argId = 1

    argId += 1
    print len(argv), argId
    if len(argv)>=argId+1: # set working folder
        workingFolder = path.abspath(argv[argId])
    else:
        workingFolder = path.abspath("./PlayGround")
    print "2\n"
    argId += 1
    if len(argv)>=argId+1: # folder to store results
        resultsFolder = path.abspath(argv[argId])
    else:
        resultsFolder = path.abspath("./RESULTS")
    print "3\n"
    argId += 1
    if len(argv)>=argId+1: # set wall time
        walltime = argv[argId]
    else:
        walltime = "%d:00:00" % (estimatedRunTime) # 1 hours per search
    print "4\n"
    argId += 1
    if len(argv)>=argId+1: # whether to compress final results folder
        compressResultsFolderAnswer = argv[argId]
    else:
        compressResultsFolderAnswer = "yes"
    print "finally: ", argId
except:
    print('Usage: generateJobs.py number_of_jobs [working_folder="./PlayGround"] [results_folder="./RESULTS"] [walltime="03:00:00" (per event)] [compress_results_folder="yes"]')
    exit()

# save config files
open("saved_configs.py", "w").writelines("""
smoothHydropackageConfigs = {
    "number_of_jobs"            :   %d,
    "working_folder"            :   "%s",
    "results_folder"            :   "%s",
    "walltime"                  :   "%s",
    "compress_results_folder"   :   "%s",
}
""" % (numberOfJobs, workingFolder, resultsFolder, walltime, compressResultsFolderAnswer)
)

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

# print welcome message
print(yellow)
greetings(3)
print(purple + "\n" + "-"*80 + "\n>>>>> Welcome to the event generator! <<<<<\n" + "-"*80 + normal)

# check prerequisites
print(green + "\n>>>>> Checking for required libraries <<<<<\n" + normal)
if not check_environment():
    print("Prerequisites not met. Install the required library first please. Aborting.")
    exit()

#check existence of executables
print(green + "\n>>>>> Checking for existence of executables <<<<<\n" + normal)
if not check_executables():
    print("Not all executables can be generated. Aborting.")
    exit()

# clean up check_prerequisites.pyc
if path.exists("check_prerequisites.pyc"): unlink("check_prerequisites.pyc")

# generate events
print(green + "\n>>>>> Generating events <<<<<\n" + normal)

# prepare directories
if not path.exists(resultsFolder): makedirs(resultsFolder)
if path.exists(workingFolder): rmtree(workingFolder)
makedirs(workingFolder)

smoothHydropackageFolder = "smoothHydropackage"
crankFolderName = "crank"
crankFolder = path.join(smoothHydropackageFolder, crankFolderName)
utilitiesFolderName = "utilities"
utilitiesFolder = utilitiesFolderName

#  parameter file into the utilities folder and backup it
copy("runHydro_shell.py", path.join(smoothHydropackageFolder, utilitiesFolder))
copy("runHydro_shell.py", resultsFolder)

# duplicate smoothHydroPackage folder to working directory, write .pbs file
for i in range(1, numberOfJobs+1):
    targetWorkingFolder = path.join(workingFolder, "node%d" % i)
    # copy folder
    copytree(smoothHydropackageFolder, targetWorkingFolder)
    open(path.join(targetWorkingFolder, "node%d.pbs" % i), "w").write(
"""
#!/usr/bin/env bash
#PBS -N smoothHydro-%d
#PBS -l walltime=%s
#PBS -l mem=8GB
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
module load python/2.7.1 # for osc oakley cluster
module load hdf5-serial
cd %s
(cd %s
    ulimit -n 1000
    python runHydro_shell.py 1> RunRecord.txt 2> ErrorRecord.txt
    cp RunRecord.txt ErrorRecord.txt ../RESULTS/
)

""" % (i, walltime, targetWorkingFolder, utilitiesFolder)
    )

print("Jobs generated. Submit them using submitJobs scripts.")



