#!/usr/bin/env python

from subprocess import call
from os import path,getcwd,makedirs,remove
from sys import stdout
import numpy as np
from glob import glob
import shutil

runHydroParameters = {
    'ecm'              :  2760,      # collision energy (GeV): 7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760
    'mode'             :  'hydro_search',   #the simulation type:  hydro[default], hybrid, hybrid_search; hydro_search, hybrid_usePrecalculated
    'model'            :  'MCGlb',   #initial condition model:  MCGlb[default], MCKLN
    'vis'              :  0.08,      #the specific shear viscosity used in the hydro simulation eta/s = 0.08 [default]
    'Tdec'             :  0.155,      #the decoupling temperature (GeV) used in the hydro simulation Tdec = 0.12 GeV [default]
    'tau0'             :  0.6,       #the hydrodynamic starting proper time (fm/c) tau0 = 0.6 fm/c [default]
    'VisBulkNorm'      :  1.0,       #the normalization factor for zeta/s(T) [default 1.0]
    'EOS'              :  's95p-v1', #s95p-v0-PCE165 [default], s95p-v1-PCE150, s95p-v1, SM-EOS-Q
    'cf_flag'          :  True,      #switch to perfrom Cooper-Frye freeze-out in pure hydro simulation cf_flag = True [default]
    'fit_flag'         :  True,      #switch to perfrom fit for normalization factor to charged multiplicity fit_flag = True [default]
    'cen'              :  '10-20',   #specify the centrality bin: All [default], e.g. 20-30
    'collision_system' :  'Pb+Pb',   #type of collision system:  Pb+Pb[default], Au+Au, Cu+Au, U+U, p+Pb, p+Au, d+Au, He+Au
    'pre_eq'           :  False,      #whether to include initial pre-equilibrium
    'parallel_mode'    :  2,         #switch to run osc2u and urqmd in parallel mode by splitting iSS resampled events to parallel_mode pieces
}

def splitParameterTable(infile_name, number_of_nodes):
    """
        split parameter table
    """
    # locate the current node
    rootDir = path.abspath("..") #assume this file is placed under utilities/
    tableLocation = path.join(rootDir, "tables")
    node_name = rootDir.split('/')[-1]  # get the name of current node
    node_index = int(node_name.split('e')[-1]) # index of tau_s for current folder
    # load original table
    params_list_data = np.loadtxt(infile_name)
    # subset data
    if params_list_data.ndim==1: # only one line of parameters
        params_selected = params_list_data 
    else:
        total_lines = params_list_data.shape[0]
        node_lines  = int(total_lines/number_of_nodes)
        params_selected = params_list_data[(node_index-1)*node_lines:(node_index*node_lines),:]
    outfile_name = path.join(rootDir, "params_list_node%d.dat"%node_index)
    np.savetxt(outfile_name, params_selected, fmt='%g',delimiter='\t')
    return params_selected


def updateParameters(params_oneline):
    """
        update parameters runHydroParameters for working in parameter search mode. 
        input: switching time, shear viscosity, switching temperature
    """
    taus, vis, tsw, visbulknorm = params_oneline
    # update dictionary
    runHydroParameters['tau0']= taus
    runHydroParameters['vis'] = vis
    runHydroParameters['Tdec']= tsw*0.001 # MeV --> GeV
    runHydroParameters['VisBulkNorm'] = visbulknorm


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

def runHydro_paramSearch():
    """
        Automatically detect the node name, extract the parameter combinations 
        for current node, and run parameter search.
    """
    number_of_nodes = len(glob(path.join('../../', 'node*')))
    params_currentNode = splitParameterTable('../tables/params_list.dat',
                                             number_of_nodes)
    if params_currentNode.ndim==1: # only one line
        params_now = params_currentNode[:]
        updateParameters(params_now)
        # form assignment string
        assignments = formAssignmentStringFromDict(runHydroParameters)
        # form executable string
        executableString = './runHydro.py' + assignments
        # execute!
        run(executableString, cwd=path.abspath("./"))
    else:        
        for i in range(params_currentNode.shape[0]):
            params_now = params_currentNode[i, :] # each line has four parameters: taus, eta/s, tdec, edec
            updateParameters(params_now)
            # form assignment string
            assignments = formAssignmentStringFromDict(runHydroParameters)
            # form executable string
            executableString = './runHydro.py' + assignments
            # execute!
            run(executableString, cwd=path.abspath("./"))


if __name__ == "__main__":
    runHydro_paramSearch()
    #runHydro_shell()