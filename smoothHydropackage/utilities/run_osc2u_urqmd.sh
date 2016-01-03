#!/usr/bin/env bash

here=$(pwd)

# run osc2u
cd $here/osc2u
./osc2u.e < OSCAR.DAT 1>run_log.dat 2>run_err.dat
rm -f OSCAR.DAT
/Users/Jia/code/smoothHydropackage/smoothHydropackage/utilities/run_osc2u_urqmd.sh
# transfer osc2u result to urqmd
mv $here/osc2u/fort.14 $here/urqmd/OSCAR.input

# run urqmd bash script
cd $here/urqmd
bash runqmd.sh 1>run_log.dat 2>run_err.dat

