#! /usr/bin/env python

import shutil
from os import path, makedirs, remove
import subprocess
from glob import glob
import numpy as np

import sys

# centrality list
cen_list = ['0-5', '5-10', '10-20', '20-30', '30-40',
            '40-50', '50-60', '60-70', '70-80']

# charged multiplicity dN/deta for 0-5% centrality
dn_deta_dict = {'5500.0': 1974.234,
                '2760.0': 1601,
                '200.0': 691,
                '62.4': 472, }
norm_factor_guess = 4.42195 # kln has flow

class color:
    """
    define colors in the terminal
    """
    purple = '\033[95m'
    cyan = '\033[96m'
    darkcyan = '\033[36m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    bold = '\033[1m'
    underline = '\033[4m'
    end = '\033[0m'

def cleanUpFolder(aDir):
    """ Delete all data files in the given directory. """
    if path.exists(aDir):
        try:
            shutil.rmtree(aDir)
        except OSError:
            pass # very likely the the folder is already empty
    makedirs(aDir)



def generate_avg_initial_condition(model, ecm, chosen_centrality, collsys,
                                   cut_type='total_entropy', pre_eq=False):
    cmd = './generateAvgprofile.py '
    args = ('-ecm %s -model %s -cen %s -cut_type %s -collision_system %s -pre_eq %s'
            % (ecm, model, chosen_centrality, cut_type, collsys, pre_eq))
    print "Generating event-averaged initial conditions..."
    print(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, cwd='./')
    p.wait()
    return

def run_pre_eq(initial_path, cen_string, run_record, err_record, tau0):
    """
        Perform pre-equilibrium evolution with averaged initial conditions
    """
    fs_path = './fs'
    fs_data_path = path.join(fs_path, 'data')
    fs_init_path = path.join(fs_data_path, 'events')
    fs_result_path = path.join(fs_data_path, 'result', 'event_1', '%g'%tau0)
    cleanUpFolder(fs_result_path)

    # prepare initial file
    shutil.copyfile('./%s/sdAvg_order_2_C%s.dat' % (initial_path, cen_string),
                    path.join(fs_init_path, 'sd_event_1_block.dat'))
    # fs
    cmd = './lm.e'
    args= (' event_mode=1 dEdyd2rdphip_dist=0 sfactor=1.0'
                + ' tau_min=%6.4f tau_max=%6.4f'
                % (tau0, tau0)) 
    sys.stdout.flush()
    run_record.write(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, stdout=run_record,
                         stderr=err_record, cwd=fs_path)
    p.wait()

def run_hydro_evo(cen_string, hydro_path, run_record, err_record,
                  norm_factor, vis, edec, tau0, pre_eq):
    """
        Perform pure hydro simulations with averaged initial conditions
    """
    initial_path = 'RESULTS/initial_conditions'
    # clean hydro initial folder
    hydro_initial_path = path.join(hydro_path, 'Initial')
    cleanUpFolder(hydro_initial_path)
    # run pre-equilibrium
    if(pre_eq == True):
        run_pre_eq(initial_path, cen_string, run_record, err_record, tau0)
        for aFile in glob(path.join('./fs/data/result/event_1/%g'%tau0, '*')):
            shutil.move(aFile, hydro_initial_path)
    else:
        shutil.copyfile('./%s/sdAvg_order_2_C%s.dat' % (initial_path, cen_string),
                    path.join(hydro_path, 'Initial', 'InitialSd.dat'))

    # hydro
    cleanUpFolder(path.join(hydro_path, 'results'))
    cmd = './VISHNew.e'
    if(pre_eq == True):
        args = (' IINIT=2 IEOS=7 iEin=0 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d'
                % (tau0, edec, vis, norm_factor, pre_eq))
    else:
        args = (' IINIT=2 IEOS=7 iEin=1 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d'
                % (tau0, edec, vis, norm_factor, pre_eq))
    print "%s : %s" % (cen_string, cmd + args)
    sys.stdout.flush()
    run_record.write(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, stdout=run_record,
                         stderr=err_record, cwd=hydro_path)
    p.wait()

def split_iSS_events(number_of_split, 
                    temp_folder = "./temp_nodes", 
                    subfolder_pattern = "node_%d",
                    temp_output_folder = "./temp_output",
                    input_file = "OSCAR.DAT"):
    """
        Sequentially implement the OSCAR events divide and conquer:
        1. split input file into N pieces using script split_events.py;
        2. replicate osc2u and UrQMD codes N times;
        3. move and rename the input file to replicated folders;
        4. run osc2u and urqmd in parallel;
        5. collect data and clean up temporary folders.
    """
    iSS_path = "./iSS"
    osc2u_path = "./osc2u"
    urqmd_path = "./urqmd"

    # split OSCAR file
    split_script = path.join(iSS_path, "split_events.py")
    # check input
    if not path.isfile(path.join(iSS_path, input_file)):
        print "No input file %s under iSS folder!"%input_file
        sys.exit(-1)
    if type(number_of_split) is not int:
        print "Check input of number_of_split: %g"%number_of_split
        number_of_split = int(number_of_split)
    split_cmd = 'python split_events.py %s %d'%(input_file, number_of_split)
    print "Start to split input file: %s"%input_file + "."*3
    p = subprocess.Popen(split_cmd, shell=True, cwd = iSS_path)
    p.wait()
    print "%s splitted!\n"%split_script

    # replicate osc2u and urqmd
    print "Start to copy osc2u, urqmd folders and input files"+'.'*3
    cleanUpFolder(temp_folder)
    for node_i in range(number_of_split):
        folder_i = path.join(temp_folder, subfolder_pattern%node_i)
        cleanUpFolder(folder_i)
        # copy codes
        shutil.copytree(osc2u_path, path.join(folder_i, 'osc2u'))
        shutil.copytree(urqmd_path, path.join(folder_i, 'urqmd'))
        # copy input file
        source_file_pattern = input_file.split('.')[0] + "_%d.dat" # according to the file format in split_events.py
        source_file = path.join(iSS_path, source_file_pattern%node_i)
        # copy osc2u and urqmd running script
        shutil.copy('./run_osc2u_urqmd.sh', folder_i)
        # check file exist
        if not path.isfile(source_file):
            print "Input file: %s does not exist! "%source_file
            sys.exit(-1)
        target_file = path.join(folder_i, "osc2u", input_file)
        shutil.move(source_file, target_file)
    print "Parallel folders created at %s!\n"%temp_folder

    # run splitted osc2u urqmd in parallel, ref to 
    process_list = []
    for node_i in range(number_of_split):
        folder_i = path.join(temp_folder, subfolder_pattern%node_i)
        # delete output files if exists
        if path.isfile(path.join(folder_i, 'osc2u', 'fort.14')):
            remove(path.join(folder_i, 'osc2u', 'fort.14'))
        if path.isfile(path.join(folder_i, 'urqmd', 'particle_list.dat')):
            remove(path.join(folder_i, 'urqmd', 'particle_list.dat'))
        cmd_string = "bash ./run_osc2u_urqmd.sh"
        p = subprocess.Popen(cmd_string, shell=True, cwd=folder_i)
        process_list.append(p)
    # waiting untill all process finish
    print "Waiting for all processes to complete, please be patient..."
    for p in process_list: p.wait()
    print "All osc2u and urqmd processes finished!"

    # collect data
    result_files = []
    cleanUpFolder(temp_output_folder)
    for node_i in range(number_of_split):
        folder_i = path.join(temp_folder, subfolder_pattern%node_i)
        source_file = path.join(folder_i, "urqmd", "particle_list.dat")
        target_file = path.join(temp_output_folder, "particle_list_%d.dat"%node_i)
        # ignore if not such file
        if not path.isfile(source_file): continue
        shutil.move(source_file, target_file)
        result_files.append(target_file)
    print "urQMD result files saved to folder %s!\n"%temp_output_folder

    # clean and quit
    shutil.rmtree(temp_folder)
    return result_files


def run_hydro_with_iS(cen_string, hydro_path, iS_path, run_record, err_record,
                      norm_factor, vis, edec, tau0, pre_eq):
    """
        Perform pure hydro simulations + Cooper Frye freeze-out
        with averaged initial conditions
    """
    run_hydro_evo(cen_string, hydro_path, run_record, err_record,
                  norm_factor, vis, edec, tau0, pre_eq)

    # move hydro results to iS
    hydro_results_path = path.join(hydro_path, 'results')
    iS_results_path = path.join(iS_path, 'results')
    cleanUpFolder(iS_results_path)
    for aFile in glob(path.join(hydro_results_path, '*')):
        if(path.exists(aFile)):
            shutil.move(aFile, iS_results_path)

    # iS
    print "%s : %s" % (cen_string, 'iS_withResonance.sh')
    sys.stdout.flush()
    p = subprocess.Popen('./iS_withResonance.sh',
                         shell=True, stdout=run_record, stderr=err_record,
                         cwd=iS_path)
    p.wait()


def run_hybrid_calculation(cen_string, model, ecm, hydro_path, iSS_path,
                           run_record, err_record,
                           norm_factor, vis, tdec, edec, tau0, eos_name,
                           pre_eq):
    """
        Perform hydro + UrQMD hybrid simulations with averaged initial
        conditions
    """
    initial_path = 'RESULTS/initial_conditions'
    result_folder = ('%s%.0fVis%gC%sTdec%gTau%g_%s'
                     % (model, ecm, vis, cen_string, tdec, tau0, eos_name))
    results_folder_path = path.join(path.abspath('./RESULTS'), result_folder)
    if path.exists(results_folder_path):
        shutil.rmtree(results_folder_path)
    makedirs(results_folder_path)

    hydro_initial_path = path.join(hydro_path, 'Initial')
    cleanUpFolder(hydro_initial_path)
    # run pre-equilibrium
    if(pre_eq == True):
        run_pre_eq(initial_path, cen_string, run_record, err_record, tau0)
        for aFile in glob(path.join('./fs/data/result/event_1/%g'%tau0, '*')):
            shutil.move(aFile, hydro_initial_path)
    else:
        shutil.copyfile('./%s/sdAvg_order_2_C%s.dat' % (initial_path, cen_string),
                    path.join(hydro_path, 'Initial', 'InitialSd.dat'))

    # hydro
    hydro_folder_path = path.join(hydro_path, 'results')
    if path.exists(path.join(hydro_folder_path)):
        shutil.rmtree(hydro_folder_path)
    makedirs(hydro_folder_path)
    cmd = './VISHNew.e'
    if(pre_eq == True):
        args = (' IINIT=2 IEOS=7 iEin=0 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d'
                % (tau0, edec, vis, norm_factor, pre_eq))
    else:
        args = (' IINIT=2 IEOS=7 iEin=1 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d'
                % (tau0, edec, vis, norm_factor, pre_eq))

    print "%s : %s" % (cen_string, cmd + args)
    sys.stdout.flush()
    run_record.write(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, stdout=run_record,
                         stderr=err_record, cwd=hydro_path)
    p.wait()

    worth_storing = []
    for aGlob in ['surface.dat', 'dec*.dat', 'ecc*.dat', 'VISH2p1_tec.dat']: #debug
        worth_storing.extend(glob(path.join(hydro_folder_path, aGlob)))
    for aFile in glob(path.join(hydro_folder_path, '*')):
        if aFile in worth_storing:
            shutil.copy(aFile, results_folder_path)

    # iSS
    iSS_folder_path = path.join(iSS_path, 'results')
    if path.exists(iSS_folder_path):
        shutil.rmtree(iSS_folder_path)
    output_file = 'OSCAR.DAT'
    if path.isfile(path.join(iSS_path, output_file)):
        remove(path.join(iSS_path, output_file))
    shutil.move(path.join(hydro_path, 'results'),
                path.join(iSS_path, 'results'))
    print "%s : %s" % (cen_string, 'iSS.e')
    sys.stdout.flush()
    p = subprocess.Popen('ulimit -n 1000; ./iSS.e', shell=True,
                         stdout=run_record, stderr=err_record, cwd=iSS_path)
    p.wait()
    worth_storing = []
    for aGlob in ['*vn*.dat']:
        worth_storing.extend(glob(path.join(iSS_folder_path, aGlob)))
    for aFile in glob(path.join(iSS_folder_path, '*')):
        if aFile in worth_storing:
            shutil.copy(aFile, results_folder_path)
    shutil.rmtree(iSS_folder_path)  # clean up

    #osc2u
    o2u_path = path.abspath('./osc2u')
    input_file = 'OSCAR.DAT'
    output_file = 'fort.14'
    if path.isfile(path.join(o2u_path, input_file)):
        remove(path.join(o2u_path, input_file))
    if path.isfile(path.join(o2u_path, output_file)):
        remove(path.join(o2u_path, output_file))
    shutil.move(path.join(iSS_path, input_file), o2u_path)
    print "%s : %s" % (cen_string, 'osu2u.e')
    sys.stdout.flush()
    p = subprocess.Popen('./osc2u.e < %s' % input_file, shell=True,
                         stdout=run_record, stderr=err_record, cwd=o2u_path)
    p.wait()
    remove(path.join(o2u_path, input_file))  # clean up

    #UrQMD
    UrQMD_path = path.abspath('./urqmd')
    input_file = 'OSCAR.input'
    output_file = 'particle_list.dat'
    if path.isfile(path.join(UrQMD_path, input_file)):
        remove(path.join(UrQMD_path, input_file))
    if path.isfile(path.join(UrQMD_path, output_file)):
        remove(path.join(UrQMD_path, output_file))
    shutil.move(path.join(o2u_path, 'fort.14'),
                path.join(UrQMD_path, input_file))
    print "%s : %s" % (cen_string, 'runqmd.sh')
    sys.stdout.flush()
    p = subprocess.Popen('bash runqmd.sh', shell=True, stdout=run_record,
                         stderr=err_record, cwd=UrQMD_path)
    p.wait()

    worth_storing = []
    for aGlob in ['particle_list.dat']:
        worth_storing.extend(glob(path.join(UrQMD_path, aGlob)))
    for aFile in glob(path.join(UrQMD_path, '*')):
        if aFile in worth_storing:
            shutil.copy(aFile, results_folder_path)
    remove(path.join(UrQMD_path, input_file))  # clean up
    remove(path.join(UrQMD_path, output_file))  # clean up


def fit_hydro(dNdeta_goal, vis, edec, tau0, pre_eq):
    """
    This function find the overall normalization factor for the hydrodynamic
    simulations at given collision energy
    :param dNdeta_goal: The measured charged hadron dN/deta in the mid rapidity
    :param vis: the specific shear viscosity for the QGP
    :param edec: the decoupling energy density in unit [GeV/fm^3]
    :param tau0: the starting time of hydrodynamic simulation in unit [fm/c]
    :return: the fitted overall normalization factor
    """
    run_record_file_name = 'run_record_fitNorm.dat'
    err_record_file_name = 'err_record_fitNorm.dat'
    run_record = open(path.join('.', run_record_file_name), 'a')
    err_record = open(path.join('.', err_record_file_name), 'a')
    hydro_path = path.abspath('./VISHNew')
    iS_path = path.abspath('./iS')
    norm_factor = norm_factor_guess
    tol = 1e-3
    target_file = 'Charged_eta_integrated_vndata.dat'
    while 1:
        icen = 0
        run_hydro_with_iS(cen_list[icen], hydro_path, iS_path, 
                          run_record, err_record,
                          norm_factor, vis, edec, tau0, pre_eq)
        # get target results
        temp_data = open(path.join(iS_path, 'results', target_file), 'r')
        dN_deta = float(temp_data.readline().split()[1])
        temp_data.close()
        print "dNdeta_goal = %g, dNdeta = %g, norm = : %g" % (
            dNdeta_goal, dN_deta, norm_factor,)
        sys.stdout.flush()
        shutil.rmtree(path.join(iS_path, 'results'))
        if abs(dN_deta - dNdeta_goal) / dNdeta_goal > tol:
            norm_factor = norm_factor * dNdeta_goal / dN_deta
        else:
            break
    shutil.move(path.join('.', run_record_file_name),
                path.abspath('./RESULTS'))
    shutil.move(path.join('.', err_record_file_name),
                path.abspath('./RESULTS'))
    return norm_factor


def run_purehydro(model, ecm, norm_factor, vis, tdec, edec, tau0,
                  eos_name, cf_flag, chosen_centrality, pre_eq):
    """
    shell function to run pure hydrodynamic simulation for all centrality bins
    """
    run_record_file_name = 'run_record_hydrowithiS.dat'
    err_record_file_name = 'err_record_hydrowithiS.dat'
    run_record = open(path.join('.', run_record_file_name), 'a')
    err_record = open(path.join('.', err_record_file_name), 'a')
    hydro_path = path.abspath('./VISHNew')
    iS_path = path.abspath('./iS')

    if chosen_centrality == 'All':
        for icen in range(len(cen_list)):
            if cf_flag:
                run_hydro_with_iS(cen_list[icen], hydro_path, iS_path,
                                  run_record, err_record,
                                  norm_factor, vis, edec, tau0, pre_eq)
                shutil.move(path.join(iS_path, 'results'),
                            path.join('RESULTS', '%s%.0fVis%gC%sTdec%gTau%g_%s'
                                      % (model, ecm, vis, cen_list[icen],
                                         tdec, tau0, eos_name)))
            else:
                run_hydro_evo(cen_list[icen], hydro_path,
                              run_record, err_record,
                              norm_factor, vis, edec, tau0, pre_eq)
                shutil.move(path.join(hydro_path, 'results'),
                            path.join('RESULTS',
                                      '%s%.0fVis%gC%sTdec%gTau%g_%s_hydroOnly'
                                      % (model, ecm, vis, cen_list[icen],
                                         tdec, tau0, eos_name)))
    else:
        if cf_flag:
            run_hydro_with_iS(chosen_centrality, hydro_path, iS_path,
                              run_record, err_record,
                              norm_factor, vis, edec, tau0, pre_eq)
            shutil.move(path.join(iS_path, 'results'),
                        path.join('RESULTS', '%s%.0fVis%gC%sTdec%gTau%g_%s'
                                  % (model, ecm, vis, chosen_centrality,
                                     tdec, tau0, eos_name)))
        else:
            run_hydro_evo(chosen_centrality, hydro_path,
                          run_record, err_record,
                          norm_factor, vis, edec, tau0, pre_eq)
            shutil.move(path.join(hydro_path, 'results'),
                        path.join('RESULTS',
                                  '%s%.0fVis%gC%sTdec%gTau%g_%s_hydroOnly'
                                  % (model, ecm, vis, chosen_centrality,
                                     tdec, tau0, eos_name)))

    shutil.move(path.join('.', run_record_file_name), 'RESULTS')
    shutil.move(path.join('.', err_record_file_name), 'RESULTS')


def run_hybrid(model, ecm, norm_factor, vis, tdec, edec,
               tau0, eos_name, chosen_centrality, pre_eq):
    """
    shell function for running hybrid calculations for all centrality bins.
    """
    run_record_file_name = 'run_record_hybrid.dat'
    err_record_file_name = 'err_record_hybrid.dat'
    run_record = open(path.join('.', run_record_file_name), 'a')
    err_record = open(path.join('.', err_record_file_name), 'a')
    hydro_path = path.abspath('./VISHNew')
    iSS_path = path.abspath('./iSS')

    if chosen_centrality == 'All':
        for icen in range(len(cen_list)):
            run_hybrid_calculation(cen_list[icen], model, ecm,
                                   hydro_path, iSS_path,
                                   run_record, err_record, norm_factor,
                                   vis, tdec, edec, tau0, eos_name, pre_eq)
    else:
        run_hybrid_calculation(chosen_centrality, model, ecm,
                               hydro_path, iSS_path,
                               run_record, err_record,
                               norm_factor, vis, tdec, edec, tau0, eos_name,
                               pre_eq)

    shutil.move(path.join('.', run_record_file_name), 'RESULTS')
    shutil.move(path.join('.', err_record_file_name), 'RESULTS')

def set_eos(eos_name, tdec):
    """
    This function replace the EOS for the whole simulation
    :param eos_name: the name of EOS
    :param tdec: the decoupling temperature (GeV)
    :return edec: the decoupling energy density (GeV/fm^3) corresponds to tdec
    """
    superMC_eos_path = './superMC/s95p-PCE'
    fs_eos_path = './fs/EOS'
    hydro_eos_path = './VISHNew/EOS/EOS_tables'
    iS_eos_path = './iS/EOS'
    iSS_eos_path = './iSS/EOS'
    if eos_name == 's95p-v0-PCE165':
        eos_files_path = './EOS/EOS_s95p/s95p_convertedtables/s95p-PCE165-v0'
    elif eos_name == 's95p-v1-PCE150':
        eos_files_path = './EOS/EOS_s95p/s95p_convertedtables/s95p-PCE-v1'
    elif eos_name == 's95p-v1':
        eos_files_path = './EOS/EOS_s95p/s95p_convertedtables/s95p-v1'
    elif eos_name == 'SM-EOS-Q':
        eos_files_path = './EOS/SMEOSQ'
    else:
        raise ValueError('invalid EOS: %s' % eos_name)

    # copy EOS files to hydro, iS, and iSS folder
    for aFile in glob(path.join(eos_files_path, '*')):
        shutil.copy(aFile, fs_eos_path)
        shutil.copy(aFile, hydro_eos_path)
        shutil.copy(aFile, iS_eos_path)
        shutil.copy(aFile, iSS_eos_path)
        
    # copy and rename EOS files to superMC folder
    shutil.copy(path.join(eos_files_path, 'EOS_PST.dat'),
                path.join(superMC_eos_path, 'EOS_converted.dat'))

    eos_file = np.loadtxt(path.join(hydro_eos_path, 'EOS_PST.dat'))
    edec = np.interp(tdec, eos_file[:,3], eos_file[:, 0])
    return edec

def run_simulations(mode, model, ecm, dN_deta, vis, tdec, tau0, eos_name,
                    cf_flag, fit_flag, chosen_centrality, collsys, pre_eq):
    """
    shell function to run simulations
    :param mode: simulation mode: hydro or hybrid
    :param model: initial condition model
    :param ecm: collision energy
    :param dN_deta: final charged multiplicity
    :param vis: the specific shear viscosity
    :param tdec: the decoupling temperature
    :param tau0: the starting time of hydrodynamic simulation
    :param eos_name: the name of EOS
    :param cf_flag: switch for Cooper-Frye freeze-out
    :param pre_eq : switch for pre-equilibrium evolution
    :return: none
    """
    print('%s mode: %s sqrt{s} = %s A GeV' % (mode, model, ecm))
    print('eta/s = %g, Tdec = %g GeV, tau0 = %g fm/c' % (vis, tdec, tau0))
    print('Pre-equilibrium: %s'%pre_eq)
    print('EOS : %s' % eos_name)

    # initial setup
    result_folder_path = './RESULTS'
    if not path.exists(result_folder_path):
        makedirs(result_folder_path)

    edec = set_eos(eos_name, tdec)

    print('preparing initial conditions ...')
    collsys_list = collsys.split('+')
    if(collsys_list[0] not in ['Au', 'Pb'] or
       collsys_list[1] not in ['Au', 'Pb']):
        modelsys = model + collsys_list[0] + collsys_list[1]
    else:
        modelsys = model

    if chosen_centrality not in cen_list and chosen_centrality != 'All':
        print("initial density profiles for %s%% centrality is not found!"
              % chosen_centrality)
        generate_flag = 'y'#raw_input("Do you want to generate one right now?")
        if generate_flag.lower() in ['yes', 'y']:
            generate_avg_initial_condition(model, ecm, chosen_centrality,
                                           collsys, pre_eq)
        else:
            sys.exit(0)
    else:
        if modelsys == "MCKLN":
            preEq_string = "_preEq" if pre_eq == True else "_nopreEq"
            initial_condition_name = '%s%.0f_sigmaNN_gauss_d0.9%s' % (modelsys, ecm, preEq_string)
        else:
            initial_condition_name = '%s%.0f_sigmaNN_gauss_d0.9' % (modelsys, ecm)            
        if path.isfile('./initial_conditions/%s.zip' % initial_condition_name):
            p = subprocess.Popen('unzip -qo %s.zip' % initial_condition_name, 
                                 shell=True, stdout=subprocess.PIPE, 
                                 cwd='./initial_conditions')
            p.wait()
            # move initial profiles
            move_cmd = "mv -f ./initial_conditions/%s/* %s/initial_conditions/"%(
                initial_condition_name, result_folder_path)
            p = subprocess.Popen(move_cmd, shell=True, stdout=subprocess.PIPE,
                cwd='./')
            p.wait()
        else:
            print("initial density profiles for %s%% centrality for %s %s " 
                  "at sqrt{s} = %g A GeV is not found!" 
                  % (chosen_centrality, model, collsys, ecm))
            generate_flag = 'y'#raw_input("Do you want to generate one right now?")
            if generate_flag.lower() in ['yes', 'y']:
                generate_avg_initial_condition(model, ecm, chosen_centrality, 
                                               collsys, pre_eq)
            else:
                sys.exit(0)

    # start to run simulations
    if fit_flag:
        print "fitting the overall normalization factor ..."
        norm_factor = fit_hydro(dN_deta, vis, edec, tau0, pre_eq)
    else:
        norm_factor = norm_factor_guess#float(input("Please input the normalization factor: "))
    if mode == 'hydro':
        print "running pure hydro simulations for centrality bin(s): %s ..."%chosen_centrality
        run_purehydro(model, ecm, norm_factor, vis, tdec, edec, tau0,
                      eos_name, cf_flag, chosen_centrality, pre_eq)
    elif mode == 'hybrid':
        print "running hybrid simulations for centrality bin(s): %s..."%chosen_centrality
        run_hybrid(model, ecm, norm_factor, vis, tdec, edec, tau0,
                   eos_name, chosen_centrality, pre_eq)
    else:
        print sys.argv[0], ': invalid running mode', mode
        sys.exit(1)


def print_help_message():
    print "Usage : "
    print(color.bold
          + "./runHydro.py -ecm ecm "
          + "[-mode mode -model model -vis vis -Tdec Tdec -tau0 tau0 "
          + "-EOS eos_name -cf_flag cf_flag -fit_flag fit_flag "
          + "-cen cen_bounds -collision_system collsys]"
          + color.end)
    print "Usage of runHydro.py command line arguments: "
    print(color.bold + "-ecm" + color.end
          + "   collision energy (GeV): "
          + color.purple + "7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760"
          + color.end)
    print(color.bold + "-mode" + color.end + "  the simulation type: "
          + color.purple + color.bold + " hydro[default]" + color.end
          + color.purple + ", hybrid" + color.end)
    print(color.bold + "-model" + color.end + " initial condition model: "
          + color.purple + color.bold + " MCGlb[default]" + color.end
          + color.purple + ", MCKLN" + color.end)
    print(color.bold + "-vis" + color.end
          + "   the specific shear viscosity used in the hydro simulation\n"
          + color.bold + "       eta/s = 0.08 [default]" + color.end)
    print(color.bold + "-Tdec" + color.end
          + "  the decoupling temperature (GeV) used in the "
          + "hydro simulation\n"
          + color.bold + "       Tdec = 0.12 GeV [default]" + color.end)
    print(color.bold + "-tau0" + color.end
          + "  the hydrodynamic starting proper time (fm/c) \n"
          + color.bold + "       tau0 = 0.6 fm/c [default]" + color.end)
    print(color.bold + "-EOS" + color.end
          + "   the equation of state for hydrodynamic simulation \n"
          + color.purple + color.bold + "       s95p-v0-PCE165 [default]"
          + color.end
          + color.purple + ", s95p-v1-PCE150, s95p-v1, SM-EOS-Q" + color.end)
    print(color.bold + "-cf_flag" + color.end
          + "   switch to perfrom Cooper-Frye freeze-out "
          + "in pure hydro simulation \n"
          + color.bold + "           cf_flag = True [default]" + color.end)
    print(color.bold + "-fit_flag" + color.end
          + "  switch to perfrom fit for normalization factor "
          + "to charged multiplicity before the simulation \n"
          + color.bold + "           fit_flag = True [default]" + color.end)
    print(color.bold + "-cen" + color.end
          + "  specify the centrality bin: "
          + color.bold + "All [default]" + color.end
          + color.purple + ', e.g. 20-30' + color.end)
    print(color.bold + "-collision_system" + color.end
          + " type of collision system: "
          + color.purple + color.bold + " Pb+Pb[default]" + color.end
          + color.purple + ", Au+Au, Cu+Au, U+U, p+Pb, p+Au, d+Au, He+Au"
          + color.end)
    print(color.bold + "-pre_eq" + color.end
          + "   switch to perfrom pre-equilibrium before hydro: "
          + color.purple + color.bold + "False [default]" + color.end)
    print(color.bold + "-h | -help" + color.end + "    This message")


if __name__ == "__main__":
    # set default values
    vis = 0.08
    tdec = 0.12  # GeV
    tau0 = 0.6  # fm/c
    mode = 'hydro'
    model = 'MCGlb'
    eos_name = 's95p-v0-PCE165'
    chosen_centrality = 'All'
    collsys = 'Pb+Pb'
    cf_flag = True
    fit_flag = True
    pre_eq = False

    while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]
        if option == '-vis':
            vis = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-Tdec':
            tdec = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-model':
            model = str(sys.argv[1])
            del sys.argv[1]
        elif option == '-tau0':
            tau0 = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-ecm':
            ecm = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-mode':
            mode = str(sys.argv[1])
            del sys.argv[1]
        elif option == '-EOS':
            eos_name = str(sys.argv[1])
            del sys.argv[1]
        elif option == '-cf_flag':
            cf_flag = (sys.argv[1] == 'True')
            del sys.argv[1]
        elif option == '-fit_flag':
            fit_flag = (sys.argv[1] == 'True')
            del sys.argv[1]
        elif option == '-cen':
            chosen_centrality = str(sys.argv[1])
            del sys.argv[1]
        elif option == '-collision_system':
            collsys = str(sys.argv[1])
            del sys.argv[1]
        elif option == '-pre_eq':
            pre_eq  = (sys.argv[1] == 'True')
            del sys.argv[1]
        elif option == '-h':
            print_help_message()
            sys.exit(0)
        else:
            print sys.argv[0], ': invalid option', option
            print_help_message()
            sys.exit(1)

    try:
        ecm_string = '%.1f' % ecm
    except NameError:
        print_help_message()
        sys.exit(1)

    # get dN/deta from the collision energy
    if 'Cu' in collsys:
        if ecm_string == '200.0':
            dN_deta = 182
        elif ecm_string == '62.4':
            dN_deta = 125.04
        else:
            print sys.argv[0], ': invalid collision energy', ecm
            sys.exit(1)
    elif 'U' in collsys:
        if ecm_string == '193.0':
            dN_deta = 1.0  # unknown yet, please use norm from Au+Au
        else:
            print sys.argv[0], ': invalid collision energy', ecm
            sys.exit(1)
    elif ecm < 62.4:
        dN_deta = 312.5 * np.log10(ecm) - 64.8
    elif ecm_string in dn_deta_dict.keys():
        dN_deta = dn_deta_dict[ecm_string]
    else:
        print sys.argv[0], ': invalid collision energy', ecm
        sys.exit(1)

    if mode in ['hydro', 'hybrid']:
        run_simulations(mode, model, ecm, dN_deta, vis, tdec, tau0, eos_name,
                        cf_flag, fit_flag, chosen_centrality, collsys, pre_eq)
    else:
        print sys.argv[0], ': invalid running mode', mode
        print_help_message()
        sys.exit(1)
