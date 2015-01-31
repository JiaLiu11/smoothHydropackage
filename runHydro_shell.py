#!/usr/bin/env python

from subprocess import call
from os import path,getcwd
from sys import stdout

runHydroParameters = {
    'ecm'              :  2760,      # collision energy (GeV): 7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760
    'mode'             :  'hydro',   #the simulation type:  hydro[default], hybrid
    'model'            :  'MCGlb',   #initial condition model:  MCGlb[default], MCKLN
    'vis'              :  0.08,      #the specific shear viscosity used in the hydro simulation eta/s = 0.08 [default]
    'Tdec'             :  0.155,      #the decoupling temperature (GeV) used in the hydro simulation Tdec = 0.12 GeV [default]
    'tau0'             :  0.6,       #the hydrodynamic starting proper time (fm/c) tau0 = 0.6 fm/c [default]
    'EOS'              :  's95p-v1', #s95p-v0-PCE165 [default], s95p-v1-PCE150, s95p-v1, SM-EOS-Q
    'cf_flag'          :  True,      #switch to perfrom Cooper-Frye freeze-out in pure hydro simulation cf_flag = True [default]
    'fit_flag'         :  True,      #switch to perfrom fit for normalization factor to charged multiplicity fit_flag = True [default]
    'cen'              :  '10-20',   #specify the centrality bin: All [default], e.g. 20-30
    'collision_system' :  'Pb+Pb',   #type of collision system:  Pb+Pb[default], Au+Au, Cu+Au, U+U, p+Pb, p+Au, d+Au, He+Au
    'pre_eq'           :  False,      #whether to include initial pre-equilibrium
}


def formAssignmentStringFromDict(aDict):
    """
        Generate a parameter-equals-value string from the given dictionary. The
        generated string has a leading blank. Extracted from iEBE package.
    """
    result = ""
    for aParameter in aDict.keys():
        result += " -{} {}".format(aParameter, aDict[aParameter])
    return result


def run(command, cwd=getcwd(), echo=True):
    """ Invoke a command from terminal and wait for it to stop. """
    if echo:
        print("-"*80)
        print("In "+cwd)
        print("Executing command: "+command)
        print("-"*80)
        stdout.flush()
    return call(command, shell=True, cwd=cwd)


def runHydro_shell():
    """

    """
    # form assignment string
    assignments = formAssignmentStringFromDict(runHydroParameters)
    # form executable string
    executableString = './runHydro.py' + assignments
    # execute!
    print executableString
    run(executableString, cwd=path.abspath("./"))


if __name__ == "__main__":
    runHydro_shell()