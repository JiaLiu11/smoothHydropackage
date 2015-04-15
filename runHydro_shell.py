#!/usr/bin/env python

from subprocess import call
from os import path,getcwd,makedirs,remove
from sys import stdout
import numpy as np
from glob import glob
import shutil

runHydroParameters = {
    'ecm'              :  2760,      # collision energy (GeV): 7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760
    'mode'             :  'hybrid_search',   #the simulation type:  hydro[default], hybrid, hybrid_search: directly fit and run in hybrid mode on chosen centrality
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


def collectObservables(flow_order_list):
    """
    Collect the particle_list to database, and only keep the analyzed database.
    Return: the results folder name
    """
    # collect to database
    results_path = path.abspath("../RESULTS")
    ebeCollector_folder = path.abspath("../EbeCollector")
    # loop over v2 and v3
    for flow_order in flow_order_list:
        params_search_log = open(path.join('..', 'param_search_log_v%d.dat'%flow_order),
            'a+')
        result_folder = ('%s%.0fVis%gC%sTdec%gTau%g_%s_v%d'
                         % (runHydroParameters["model"], runHydroParameters["ecm"], 
                            runHydroParameters["vis"], runHydroParameters["cen"], 
                            runHydroParameters["Tdec"], runHydroParameters["tau0"], 
                            runHydroParameters["EOS"], flow_order))
        result_path_now = path.join(results_path, result_folder)
        particle_list_files = glob(path.join(result_path_now, 'particle_list*.dat'))
        # decide if parallel model ends properly
        if runHydroParameters["parallel_mode"]!=len(particle_list_files):
            print "Warning: parallel run mode does not end properly!"
            print "only %d of %d output file exist!\n"%(len(particle_list_files),
                                                runHydroParameters["parallel_mode"])
        # move files to EbeCollector folder
        if path.exists(path.join(ebeCollector_folder, 'event-1')):
            shutil.rmtree(path.join(ebeCollector_folder, 'event-1'))
            makedirs(path.join(ebeCollector_folder, 'event-1'))
        for aFile in particle_list_files:
            shutil.move(aFile, path.join(ebeCollector_folder, 'event-1'))
        # start to collect database --> from particle_list to particles.db
        if len(particle_list_files)==1:
            call("python EbeCollectorShell_particlesUrQMD.py ./", 
                shell=True, cwd=ebeCollector_folder)
        else:
            call("python EbeCollectorShell_particlesUrQMD_parallel.py ./ %d"%len(particle_list_files),
                shell=True, cwd=ebeCollector_folder)
        # silt particle info  --> from particles.db to analyzed_particles.db
        call("python particleReaderShell.py particles.db",
            shell=True, cwd=ebeCollector_folder)
        # output final observables
        call("python AnalyzedEventsReader.py analyzed_particles.db",
            shell=True, cwd=ebeCollector_folder)
        # collect data
        params_output = np.loadtxt(path.join(ebeCollector_folder, 
                                            'paramSearch_result.dat'))
        params_search_log.write(" ".join(map(lambda(x): '%10.8e'%x, params_output[:]))+'\n')
        params_search_log.close()
        # save the analyzed database to result folder
        shutil.move(path.join(ebeCollector_folder, 'analyzed_particles.db'),
            result_path_now)
        # clean up source files
        for aFile in particle_list_files:
            if path.isfile(path.join(ebeCollector_folder, 'event-1', aFile)):
                remove(path.join(ebeCollector_folder, 'event-1', aFile))
        remove(path.join(ebeCollector_folder, 'particles.db'))
    return result_folder


def updateParameters(params_oneline):
    """
        update parameters runHydroParameters for working in parameter search mode. 
        input: switching time, shear viscosity, switching temperature
    """
    taus, vis, tsw = params_oneline
    # update dictionary
    runHydroParameters['tau0']= taus
    runHydroParameters['vis'] = vis
    runHydroParameters['Tdec']= tsw*0.001 # MeV --> GeV


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
    number_of_nodes = len(glob(path.join('../../', 'node?')))
    params_currentNode = splitParameterTable('../tables/params_list.dat',
                                             number_of_nodes)
    backup_path = path.join('/nfs/gpfs/PAS0254/paramSearch',
        '%s_%d'%(runHydroParameters['model'], runHydroParameters['pre_eq']))
    # if not path.exists(backup_path):
    #     makedirs(backup_path)
    flow_order_list = [2,3]
    if params_currentNode.ndim==1: # only one line
        params_now = params_currentNode[:-1]
        updateParameters(params_now)
        # form assignment string
        assignments = formAssignmentStringFromDict(runHydroParameters)
        # form executable string
        executableString = './runHydro.py' + assignments
        # execute!
        run(executableString, cwd=path.abspath("./"))
        # collect
        save_folder_name = collectObservables(flow_order_list)
        # compress
        zipped_file_name = save_folder_name.split('_v')[0]+'.zip'
        zip_cmd = ('zip -r -q -m %s'%zipped_file_name + 
            ' ./MC*')
        print "runHydro_shell: start to compress file: %s......"%zipped_file_name
        call(zip_cmd, shell=True, cwd=path.abspath('../RESULTS'))
        # backup 
        if path.exists(path.join(backup_path, zipped_file_name)):
            remove(path.exists(backup_path, zipped_file_name))
        shutil.move(path.join('../RESULTS', zipped_file_name),
            backup_path)
        print "runHydro_shell: file %s saved to project folder!\n"%zipped_file_name
    else:        
        for i in range(params_currentNode.shape[0]):
            params_now = params_currentNode[i, :-1] # each line has four parameters: taus, eta/s, tdec, edec
            updateParameters(params_now)
            # form assignment string
            assignments = formAssignmentStringFromDict(runHydroParameters)
            # form executable string
            executableString = './runHydro.py' + assignments
            # execute!
            run(executableString, cwd=path.abspath("./"))
            # collect
            save_folder_name = collectObservables(flow_order_list)
            # compress
            zipped_file_name = save_folder_name.split('_v')[0]+'.zip'
            zip_cmd = ('zip -r -q -m %s'%zipped_file_name + 
                ' ./MC*')
            print "runHydro_shell: start to compress file: %s......"%zipped_file_name
            call(zip_cmd, shell=True, cwd=path.abspath('../RESULTS'))
            # backup 
            if path.exists(path.join(backup_path, zipped_file_name)):
                remove(path.exists(backup_path, zipped_file_name))
            shutil.move(path.join('../RESULTS', zipped_file_name),
                backup_path)
            print "runHydro_shell: file %s saved to project folder!\n"%zipped_file_name


if __name__ == "__main__":
    runHydro_paramSearch()
    #runHydro_shell()