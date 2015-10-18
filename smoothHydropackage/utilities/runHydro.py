#! /usr/bin/env python

import shutil
from os import path, makedirs, remove
import subprocess
from glob import glob
import numpy as np
import re
import sys
import time

# centrality list
cen_list = ['0-5', '5-10', '10-20', '20-30', '30-40',
            '40-50', '50-60', '60-70', '70-80']
dn_deta_ECM2760_dict = {
  '0-5'         :        1601, 
  '5-10'        :        1294,
  '10-20'       :        966,
  '20-30'       :        649,
  '30-40'       :        426,
  '40-50'       :        261,
  '50-60'       :        149,
  '60-70'       :        76,
  '70-80'       :        35,
}

# charged multiplicity dN/deta for 0-5% centrality
dn_deta_dict = {'5500.0': 1974.234,
                '2760.0': 1601,
                '200.0': 691,
                '62.4': 472, }
rootDir = path.abspath('../')
project_directory = '/nfs/gpfs/PAS0254/paramSearch'

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

def linearFitSfactor(tau_s, eta_s, t_sw, model="MCGlb", pre_eq=False, visbulknorm=0.0):
    """
        Use the linear fit result to give an guess of the scaling factor for given
        switching time, shear viscosity and switching temperature.
        The linear fit has the form:
        sfactor = A0 + A1*tau_s + A2*eta_s + A3*t_sw + A4*tau_s^2 + A5*tau_s_3 
                   + A6*eta_s^2 + A7*tau_s*eta_s
    """
    model_str = model+"_%d"%pre_eq
    sfactor_lm = -1
    if abs(visbulknorm)<1e-6:
        # choose linear coefficients
        if model_str == "MCGlb_0":
            coeff_list = np.array([192.14190497,-284.36945039,-201.68632700,-0.04096829,176.83509284,
                -39.69363878,-441.47816092,178.69094531])
        elif model_str == "MCGlb_1":
            coeff_list = np.array([36.02240022,-7.32523631,-138.53732379,-0.00868909,-3.84600983,
                1.65251890,174.27023844,67.81849622])
        elif model_str == "MCKLN_0":
            coeff_list = np.array([44.410314428,-66.071626069,-92.632445045,-0.003807115,46.914963646,
                -14.722666889,48.305633428,67.176481583])
        elif model_str == "MCKLN_1":
            coeff_list = np.array([14.16295555,3.90769983,-36.13522848,-0.00381599,-9.33504831,
                3.46467929,28.07882091,17.86500695])
        else:
            print "linearFitSfactor: cannot predict scaling factor for run mode: %s, pre-eq. = %s"%(model, pre_eq)
            return norm_factor_default # default value
        # change unit: GeV --> MeV
        t_sw = t_sw*1000.0
        # combine to linear model
        sfactor_lm = (coeff_list[0] + coeff_list[1]*tau_s + coeff_list[2]*eta_s + coeff_list[3]*t_sw
                + coeff_list[4]*tau_s**2.0 + coeff_list[5]*tau_s**3.0 
                + coeff_list[6]*eta_s**2.0 + coeff_list[7]*tau_s*eta_s)
    else:
        # choose linear coefficients
        if model_str == "MCGlb_0":
            coeff_list = np.array([1,0,0,0,0,0,0,0])
        elif model_str == "MCGlb_1":
            coeff_list = np.array([32.4855955,-15.8834285,-54.5494184,-3.7276310,7.9139809,-1.6334713,
                29.3104912,0.6376767])
        elif model_str == "MCKLN_0":
            coeff_list = np.array([35.291801, -52.281676,-33.673811,-2.321871,34.195697,-8.403409,
                21.454700, 1.162961])
        elif model_str == "MCKLN_1":
            coeff_list = np.array([12.5669,-2.8823,-10.8335,-1.5868,0.80698,-0.066272,4.83818,
                0.1713826])
        # combine to linear model
        sfactor_lm = (coeff_list[0] + coeff_list[1]*tau_s + coeff_list[2]*eta_s + coeff_list[3]*visbulknorm
                + coeff_list[4]*tau_s**2.0 + coeff_list[5]*tau_s**3.0 
                + coeff_list[6]*tau_s*eta_s + coeff_list[7]*tau_s*visbulknorm)
    return sfactor_lm

def generate_avg_initial_condition(model, ecm, chosen_centrality, collsys,
                                   cut_type='total_entropy', pre_eq=False):
    cmd = './generateAvgprofile.py '
    args = ('-ecm %s -model %s -cen %s -cut_type %s -collision_system %s -pre_eq %s'
            % (ecm, model, chosen_centrality, cut_type, collsys, pre_eq))
    print "Generating event-averaged initial conditions..."
    print(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, cwd=rootDir)
    p.wait()
    return

def run_pre_eq(initial_path, cen_string, run_record, err_record, tau0, flow_order=2):
    """
        Perform pre-equilibrium evolution with averaged initial conditions
    """
    fs_path = path.join(rootDir, 'fs')
    fs_data_path = path.join(fs_path, 'data')
    fs_init_path = path.join(fs_data_path, 'events')
    fs_result_path = path.join(fs_data_path, 'result', 'event_1', '%g'%tau0)
    cleanUpFolder(fs_result_path)

    # prepare initial file
    if not flow_order==2:
        shutil.copyfile('%s/sdAvg_order_%d_C%s.dat' % (initial_path,flow_order, cen_string),
                        path.join(fs_init_path, 'sd_event_1_block.dat'))
    else:
        shutil.copyfile('%s/sdAvg_order_2_C%s.dat' % (initial_path, cen_string),
                        path.join(fs_init_path, 'sd_event_1_block.dat'))
    # fs
    cmd = './lm.e'
    args= (' event_mode=1 dEdyd2rdphip_dist=0 sfactor=1.0'
                + ' taumin=%g taumax=%g'
                % (tau0, tau0)) 
    sys.stdout.flush()
    run_record.write(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, stdout=run_record,
                         stderr=err_record, cwd=fs_path)
    p.wait()

def run_hydro_evo(cen_string, hydro_path, run_record, err_record,
                  norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq):
    """
        Perform pure hydro simulations with averaged initial conditions
    """
    initial_path = path.join(rootDir, 'RESULTS/initial_conditions')
    # clean hydro initial folder
    hydro_initial_path = path.join(hydro_path, 'Initial')
    cleanUpFolder(hydro_initial_path)
    # run pre-equilibrium
    if(pre_eq == True):
        run_pre_eq(initial_path, cen_string, run_record, err_record, tau0)
        for aFile in glob(path.join(rootDir, 'fs/data/result/event_1/%g'%tau0, '*')):
            shutil.move(aFile, hydro_initial_path)
    else:
        shutil.copyfile('%s/sdAvg_order_2_C%s.dat' % (initial_path, cen_string),
                    path.join(hydro_path, 'Initial', 'InitialSd.dat'))

    # hydro
    cleanUpFolder(path.join(hydro_path, 'results'))
    cmd = './VISHNew.e'
    vis_temp = vis
    if vis<1e-3: vis_temp = 0
    if(pre_eq == True):
        args = (' IINIT=2 IEOS=7 iEin=0 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d visbulknorm=%.6f'
                % (tau0, edec, vis_temp, norm_factor, pre_eq, VisBulkNorm))
    else:
        args = (' IINIT=2 IEOS=7 iEin=1 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d visbulknorm=%.6f'
                % (tau0, edec, vis_temp, norm_factor, pre_eq, VisBulkNorm))
    print "%s : %s" % (cen_string, cmd + args)
    sys.stdout.flush()
    run_record.write(cmd + args)
    p = subprocess.Popen(cmd + args, shell=True, stdout=run_record,
                         stderr=err_record, cwd=hydro_path)
    p.wait()

def split_iSS_events(number_of_split, 
                    temp_folder = path.join(rootDir, "temp_nodes"), 
                    subfolder_pattern = "node_%d",
                    output_folder = path.join(rootDir, "temp_output"),
                    input_file = "OSCAR.DAT"):
    """
        Sequentially implement the OSCAR events divide and conquer:
        1. split input file into N pieces using script split_events.py;
        2. replicate osc2u and UrQMD codes N times;
        3. move and rename the input file to replicated folders;
        4. run osc2u and urqmd in parallel;
        5. collect data and clean up temporary folders.
    """
    iSS_path = path.join(rootDir, "iSS")
    osc2u_path = path.join(rootDir, "osc2u")
    urqmd_path = path.join(rootDir, "urqmd")

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
    remove(path.join(iSS_path, input_file)) # delete OSCAR.DAT to save disk space

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
        shutil.copy('run_osc2u_urqmd.sh', folder_i)
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
    urqmd_start = time.time()
    for p in process_list: p.wait()
    print "All osc2u and urqmd processes finished!"
    urqmd_end = time.time()
    print "UrQMD running time: %f seconds!"%(urqmd_end-urqmd_start)

    # collect data
    result_files = []
    if not path.exists(output_folder): makedirs(output_folder)
    for node_i in range(number_of_split):
        folder_i = path.join(temp_folder, subfolder_pattern%node_i)
        source_file = path.join(folder_i, "urqmd", "particle_list.dat")
        target_file = path.join(output_folder, "particle_list_%d.dat"%node_i)
        # ignore if not such file
        if not path.isfile(source_file): continue
        # overwrite if target file exists
        if path.isfile(target_file): remove(target_file)
        shutil.move(source_file, target_file)
        result_files.append(target_file)
    print "urQMD result files saved to folder %s!\n"%output_folder

    # clean and quit
    shutil.rmtree(temp_folder)
    return result_files


def collectObservables(result_folder, parallel_mode):
    """
    Collect the particle_list to database, and only keep the analyzed database.
    Return: the results folder name
    """
    # extract run parameters for result_folder
    num_from_str = map(float, re.findall(r"[-+]?\d*\.\d+|\d+",result_folder))
    taus, etas, tdec, VisBulkNorm = [num_from_str[i] for i in [5,1,4,6]]

    # collect to database
    results_path = path.join(rootDir, "RESULTS")
    ebeCollector_folder = path.join(rootDir,"EbeCollector")
    results_folder_path = path.join(results_path, result_folder)
    flow_order = result_folder.split('_v')[-1]
    params_search_log = open(path.join('..', 'param_search_log_v%s.dat'%flow_order),
        'a+')
    particle_list_files = glob(path.join(results_folder_path, 'particle_list*.dat'))
    # decide if parallel model ends properly
    if parallel_mode!=len(particle_list_files):
        print "Warning: parallel run mode does not end properly!"
        print "only %d of %d output file exist!\n"%(len(particle_list_files),
                                            parallel_mode)
    # move files to EbeCollector folder
    if path.exists(path.join(ebeCollector_folder, 'event-1')):
        shutil.rmtree(path.join(ebeCollector_folder, 'event-1'))
    makedirs(path.join(ebeCollector_folder, 'event-1'))
    print "collectMoveDB"+"="*20
    for aFile in particle_list_files:
        shutil.move(aFile, path.join(ebeCollector_folder, 'event-1'))
    # start to collect database --> from particle_list to particles.db
    if len(particle_list_files) == 1:
        subprocess.call("python EbeCollectorShell_particlesUrQMD.py ./ 1>run_log.dat 2>run_err.dat", 
            shell=True, cwd=ebeCollector_folder)
    else:
        subprocess.call("python EbeCollectorShell_particlesUrQMD_parallel.py ./ %d 1>run_log.dat 2>run_err.dat"%len(particle_list_files),
            shell=True, cwd=ebeCollector_folder)
    # silt particle info  --> from particles.db to analyzed_particles.db
    subprocess.call("python particleReaderShell.py particles.db 1>run_log.dat 2>run_err.dat",
        shell=True, cwd=ebeCollector_folder)
    # output final observables
    subprocess.call("python AnalyzedEventsReader.py analyzed_particles.db 1>run_log.dat 2>run_err.dat",
        shell=True, cwd=ebeCollector_folder)
    # collect data
    params_output = np.loadtxt(path.join(ebeCollector_folder, 
                                        'paramSearch_result.dat'))
    result_line = ("%.8f \t %.8f \t %.8f \t %.8f \t"%(taus, etas, tdec, VisBulkNorm)+
        " ".join(map(lambda(x): '%10.8e'%x, params_output[:]))+'\n')
    params_search_log.write(result_line)
    params_search_log.close()
    # save the analyzed database to result folder
    shutil.move(path.join(ebeCollector_folder, 'analyzed_particles.db'),
        results_folder_path)
    # shutil.move(path.join(ebeCollector_folder, 'particles.db'),
    #     results_folder_path)
    # clean up source files
    for aFile in particle_list_files:
        if path.isfile(path.join(ebeCollector_folder, 'event-1', aFile)):
            remove(path.join(ebeCollector_folder, 'event-1', aFile))
    return results_folder_path


def moveDB(source_folder, model, pre_eq):
    backup_path = path.join(project_directory, 
        '%s_%d'%(model, pre_eq))
    # compress
    zipped_file_name = source_folder+'.zip' 
    zip_cmd = ('zip -r -q -m %s'%zipped_file_name + 
        ' %s/'%source_folder)
    print "Start to compress file: %s......"%zipped_file_name
    subprocess.call(zip_cmd, shell=True, cwd=path.join(rootDir, 'RESULTS'))
    if not path.exists(backup_path):
        print "Cannot find backup path: %s"%backup_path
    else:
        # backup 
        if path.exists(path.join(backup_path, zipped_file_name)):
            remove(path.join(backup_path, zipped_file_name))
        shutil.move(path.join(rootDir, 'RESULTS', zipped_file_name),
            backup_path)
        print "File %s saved to project folder!"%zipped_file_name
    print "collectMoveDB"+"="*20+'\n'
    if path.isfile(path.join(backup_path, zipped_file_name)):
        print "parameter search file saved!"
        flag = True
    else:
        flag = False
    return flag

def backupiSSOSCAR(iSS_location, backup_file_name):
    backup_path = path.join(project_directory,
        '%s_%d'%(model, pre_eq), 'iSS_backup')
    zipped_file_name = '%s.zip'%backup_file_name
    # compress
    zip_cmd = 'zip -q %s.zip OSCAR.DAT'%backup_file_name
    print "Start to compress iSS file: %s......"%backup_file_name
    subprocess.call(zip_cmd, shell=True, cwd=iSS_location)
    if not path.exists(backup_path):
        print "Cannot find backup path: %s, creating it..."%backup_path
        makedirs(backup_path)
    # backup 
    if path.exists(path.join(backup_path, zipped_file_name)):
        remove(path.join(backup_path, zipped_file_name))
    shutil.move(path.join(iSS_location, zipped_file_name),
        backup_path)
    print "File %s saved to project iSS folder!"%zipped_file_name


def run_hydro_with_iS(cen_string, hydro_path, iS_path, run_record, err_record,
                      norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq):
    """
        Perform pure hydro simulations + Cooper Frye freeze-out
        with averaged initial conditions
    """
    run_hydro_evo(cen_string, hydro_path, run_record, err_record,
                  norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq)

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
                           norm_factor, vis, tdec, edec, tau0, VisBulkNorm, eos_name,
                           pre_eq, parallel_mode, flow_order=2):
    """
        Perform hydro + UrQMD hybrid simulations with averaged initial
        conditions
    """
    initial_path = path.join(rootDir, 'RESULTS/initial_conditions')
    if not flow_order==2:
        result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_v%d'
                         % (model, ecm, vis, cen_string, tdec, tau0, VisBulkNorm, 
                            eos_name, flow_order))
    else:
        result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s'
                         % (model, ecm, vis, cen_string, tdec, tau0, VisBulkNorm,
                            eos_name))
    results_folder_path = path.join(rootDir, 'RESULTS', result_folder)
    if path.exists(results_folder_path):
        shutil.rmtree(results_folder_path)
    makedirs(results_folder_path)

    hydro_initial_path = path.join(hydro_path, 'Initial')
    cleanUpFolder(hydro_initial_path)
    # run pre-equilibrium
    if(pre_eq == True):
        run_pre_eq(initial_path, cen_string, run_record, err_record, tau0, flow_order)
        for aFile in glob(path.join(rootDir, 'fs/data/result/event_1/%g'%tau0, '*')):
            shutil.move(aFile, hydro_initial_path)
    else:
        if not flow_order==2:
            shutil.copyfile('%s/sdAvg_order_%d_C%s.dat' % (initial_path, flow_order, cen_string),
                    path.join(hydro_path, 'Initial', 'InitialSd.dat'))
        else:
            shutil.copyfile('%s/sdAvg_order_2_C%s.dat' % (initial_path, cen_string),
                        path.join(hydro_path, 'Initial', 'InitialSd.dat'))

    # hydro
    hydro_folder_path = path.join(hydro_path, 'results')
    if path.exists(path.join(hydro_folder_path)):
        shutil.rmtree(hydro_folder_path)
    makedirs(hydro_folder_path)
    cmd = './VISHNew.e'
    if(pre_eq == True):
        args = (' IINIT=2 IEOS=7 iEin=0 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d visbulknorm=%.6f'
                % (tau0, edec, vis, norm_factor, pre_eq, VisBulkNorm))
    else:
        args = (' IINIT=2 IEOS=7 iEin=1 iLS=200'
                + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d visbulknorm=%.6f'
                % (tau0, edec, vis, norm_factor, pre_eq, VisBulkNorm))

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

    if parallel_mode > 1:
        # run subsequent programs in parallel
        result_files = split_iSS_events(number_of_split = parallel_mode,
                                        output_folder = results_folder_path)
    else:
        #osc2u
        o2u_path = path.join(rootDir,'osc2u')
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
        UrQMD_path = path.abspath(rootDir,'urqmd')
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


def fit_hydro(dNdeta_goal, vis, edec, tau0, VisBulkNorm, pre_eq, norm_factor_guess=10.0, cen_string="0-5"):
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
    run_record = open(path.join(rootDir, run_record_file_name), 'w+')
    err_record = open(path.join(rootDir, err_record_file_name), 'w+')
    hydro_path = path.join(rootDir, 'VISHNew')
    iS_path = path.join(rootDir, 'iS')
    norm_factor = norm_factor_guess
    tol = 1e-2 # approximately 10 when dNdeta_goal = 1600
    target_file = 'Charged_eta_integrated_vndata.dat'
    sfactor_log = open(path.join('..', 'sfactor_log.dat'), 'a+')
    while 1:
        try:
            icen = cen_list.index(cen_string)
            print "fit_hydro: fitting on centrality %s%%"%cen_string
        except:
            print "No such centrality: %s%%"%cen_string
            sys.exit(0)

        run_hydro_with_iS(cen_list[icen], hydro_path, iS_path, 
                          run_record, err_record,
                          norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq)
        # get target results
        temp_data = open(path.join(iS_path, 'results', target_file), 'r')
        dN_deta = float(temp_data.readline().split()[1])
        temp_data.close()
        print "dNdeta_goal = %g, dNdeta = %g, norm = : %g" % (
            dNdeta_goal, dN_deta, norm_factor,)
        sys.stdout.flush()
        if abs(dN_deta - dNdeta_goal) / dNdeta_goal > tol:
            sfactor_log.write("0 %.6f \t %.6f \t %.6f\t  %.6f  \t  %10.6e \t   %10.8f  \n"%(
                tau0, vis, edec, VisBulkNorm, norm_factor, dN_deta))
            norm_factor = norm_factor * dNdeta_goal / dN_deta
            shutil.rmtree(path.join(iS_path, 'results')) 
        else:
            sfactor_log.write("1 %.6f \t %.6f \t %.6f\t  %.6f  \t  %10.6e \t   %10.8f  \n"%(
                tau0, vis, edec, VisBulkNorm, norm_factor, dN_deta))
            break
    sfactor_log.close()
    if(path.isfile(path.join(rootDir, 'RESULTS', run_record_file_name))):
        remove(path.join(rootDir, 'RESULTS', run_record_file_name))
    if(path.isfile(path.join(rootDir, 'RESULTS', err_record_file_name))):
        remove(path.join(rootDir, 'RESULTS', err_record_file_name))    
    shutil.move(path.join(rootDir, run_record_file_name),
                path.join(rootDir, 'RESULTS'))
    shutil.move(path.join(rootDir, err_record_file_name),
                path.join(rootDir, 'RESULTS'))
    return norm_factor


def run_purehydro(model, ecm, norm_factor, vis, tdec, edec, tau0, VisBulkNorm,
                  eos_name, cf_flag, chosen_centrality, pre_eq):
    """
    shell function to run pure hydrodynamic simulation for all centrality bins
    """
    run_record_file_name = 'run_record_hydrowithiS.dat'
    err_record_file_name = 'err_record_hydrowithiS.dat'
    run_record = open(path.join(rootDir, run_record_file_name), 'a')
    err_record = open(path.join(rootDir, err_record_file_name), 'a')
    hydro_path = path.join(rootDir, 'VISHNew')
    iS_path = path.join(rootDir, 'iS')

    if chosen_centrality == 'All':
        for icen in range(len(cen_list)):
            if cf_flag:
                run_hydro_with_iS(cen_list[icen], hydro_path, iS_path,
                                  run_record, err_record,
                                  norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq)
                shutil.move(path.join(iS_path, 'results'),
                            path.join(rootDir,'RESULTS', '%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s'
                                      % (model, ecm, vis, cen_list[icen],
                                         tdec, tau0, VisBulkNorm, eos_name)))
            else:
                run_hydro_evo(cen_list[icen], hydro_path,
                              run_record, err_record,
                              norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq)
                shutil.move(path.join(hydro_path, 'results'),
                            path.join(rootDir, 'RESULTS',
                                      '%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_hydroOnly'
                                      % (model, ecm, vis, cen_list[icen],
                                         tdec, tau0, VisBulkNorm, eos_name)))
    else:
        if cf_flag:
            run_hydro_with_iS(chosen_centrality, hydro_path, iS_path,
                              run_record, err_record,
                              norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq)
            shutil.move(path.join(iS_path, 'results'),
                        path.join(rootDir, 'RESULTS', '%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s'
                                  % (model, ecm, vis, chosen_centrality,
                                     tdec, tau0, VisBulkNorm, eos_name)))
        else:
            run_hydro_evo(chosen_centrality, hydro_path,
                          run_record, err_record,
                          norm_factor, vis, edec, tau0, VisBulkNorm, pre_eq)
            shutil.move(path.join(hydro_path, 'results'),
                        path.join(rootDir, 'RESULTS',
                                  '%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_hydroOnly'
                                  % (model, ecm, vis, chosen_centrality,
                                     tdec, tau0, VisBulkNorm, eos_name)))
            
    if(path.isfile(path.join(rootDir, 'RESULTS', run_record_file_name))):
        remove(path.join(rootDir, 'RESULTS', run_record_file_name))
    if(path.isfile(path.join(rootDir, 'RESULTS', err_record_file_name))):
        remove(path.join(rootDir, 'RESULTS', err_record_file_name))  
    shutil.move(path.join(rootDir, run_record_file_name), 
        path.join(rootDir, 'RESULTS'))
    shutil.move(path.join(rootDir, err_record_file_name), 
        path.join(rootDir, 'RESULTS'))


def run_hybrid(model, ecm, norm_factor, vis, tdec, edec,
               tau0, VisBulkNorm, eos_name, chosen_centrality, pre_eq,
               parallel_mode):
    """
    shell function for running hybrid calculations for all centrality bins.
    """
    run_record_file_name = 'run_record_hybrid.dat'
    err_record_file_name = 'err_record_hybrid.dat'
    run_record = open(path.join(rootDir, run_record_file_name), 'a')
    err_record = open(path.join(rootDir, err_record_file_name), 'a')
    hydro_path = path.join(rootDir, 'VISHNew')
    iSS_path = path.join(rootDir, 'iSS')

    if chosen_centrality == 'All':
        for icen in range(len(cen_list)):
            run_hybrid_calculation(cen_list[icen], model, ecm,
                                   hydro_path, iSS_path,
                                   run_record, err_record, norm_factor,
                                   vis, tdec, edec, tau0, VisBulkNorm, eos_name, pre_eq,
                                   parallel_mode)
    else:
        run_hybrid_calculation(chosen_centrality, model, ecm,
                               hydro_path, iSS_path,
                               run_record, err_record,
                               norm_factor, vis, tdec, edec, tau0, VisBulkNorm, eos_name,
                               pre_eq, parallel_mode)

    shutil.copy(path.join(rootDir, run_record_file_name), 
        path.join(rootDir, 'RESULTS'))
    shutil.copy(path.join(rootDir, err_record_file_name), 
        path.join(rootDir, 'RESULTS'))

def run_hybrid_search(model, ecm, norm_factor, vis, tdec, edec,
               tau0, VisBulkNorm, eos_name, chosen_centrality, pre_eq,
               parallel_mode):
    """
    Shell for parameter search mode simulation. Divert iS input to iSS 
    and start urqmd thereafter.
    """
    run_record_file_name = 'run_record_hybrid_search.dat'
    err_record_file_name = 'err_record_hybrid_search.dat'
    run_record = open(path.join(rootDir, run_record_file_name), 'w+')
    err_record = open(path.join(rootDir, err_record_file_name), 'w+')
    hydro_path = path.join(rootDir, 'VISHNew')
    iS_results_path  = path.join(rootDir, 'iS', 'results')
    iSS_path = path.join(rootDir, 'iSS')

    # run after burner for v2
    result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_v%d'
                     % (model, ecm, vis, chosen_centrality, tdec, tau0, VisBulkNorm,
                        eos_name, 2))
    results_folder_path = path.join(rootDir, 'RESULTS', result_folder)
    if path.exists(results_folder_path):
        shutil.rmtree(results_folder_path)
    makedirs(results_folder_path)

    # save files after hydro run
    worth_storing = []
    for aGlob in ['surface.dat', 'dec*.dat', 'ecc*.dat', 'VISH2p1_tec.dat']: #debug
        worth_storing.extend(glob(path.join(iS_results_path, aGlob)))
    for aFile in glob(path.join(iS_results_path, '*')):
        if aFile in worth_storing:
            shutil.copy(aFile, results_folder_path)
    run_afterBurner(iS_results_path, chosen_centrality, run_record, err_record,
        results_folder_path, parallel_mode)
    # collect db and backup v2 search result
    collectObservables(result_folder, parallel_mode)
    moveDB(result_folder, model, pre_eq)

    # run the search for v3
    print "Start to search v3:"
    run_hybrid_calculation(chosen_centrality, model, ecm,
                           hydro_path, iSS_path,
                           run_record, err_record,
                           norm_factor, vis, tdec, edec, tau0, VisBulkNorm, eos_name,
                           pre_eq, parallel_mode, 3)
    result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_v%d'
                     % (model, ecm, vis, chosen_centrality, tdec, tau0, VisBulkNorm, 
                        eos_name, 3))
    # collect db and backup v3 search result
    collectObservables(result_folder, parallel_mode)
    moveDB(result_folder, model, pre_eq)

    run_record.close()
    err_record.close()
    shutil.copy(path.join(rootDir, run_record_file_name), 
        path.join(results_folder_path))
    shutil.copy(path.join(rootDir, err_record_file_name), 
        path.join(results_folder_path))


def run_hydro_search(model, ecm, norm_factor, vis, tdec, edec,
               tau0, VisBulkNorm, eos_name, chosen_centrality, pre_eq):
    """
    Shell for parameter search mode simulation in pure hydro. 
    Save v2 data from fit_hydro, and run v3 in pure hydro mode without iS.
    """
    initial_path = path.join(rootDir, 'RESULTS/initial_conditions')
    hydro_path = path.join(rootDir, 'VISHNew')
    iS_path = path.join(rootDir, 'iS')
    iS_results_path  = path.join(rootDir, 'iS', 'results')
    iSS_path = path.join(rootDir, 'iSS')

    run_record_file_name = 'run_record_hydro_search.dat'
    err_record_file_name = 'err_record_hydro_search.dat'
    run_record = open(path.join(rootDir, run_record_file_name), 'w+')
    err_record = open(path.join(rootDir, err_record_file_name), 'w+')

    # prepare data backup folder for v2
    flow_order = 2
    result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_v%d'
                     % (model, ecm, vis, chosen_centrality, tdec, tau0, VisBulkNorm,
                        eos_name, flow_order))
    if path.isfile(path.join(project_directory, 
                            '%s_%d'%(model, pre_eq),
                            'iSS_backup', '%s.zip'%result_folder)):
        print "hydro search results have already been found, skipping v2..."
    else:
        print "start to generate iSS events for v2 ... "
        results_folder_path = path.join(rootDir, 'RESULTS', result_folder)
        if path.exists(results_folder_path):
            shutil.rmtree(results_folder_path)
        makedirs(results_folder_path)
        # save v2 data - only hydro output
        worth_storing = []
        for aGlob in ['surface.dat', 'dec*.dat', 'ecc*.dat', 'VISH2p1_tec.dat']:
            worth_storing.extend(glob(path.join(iS_results_path, aGlob)))
        for aFile in glob(path.join(iS_results_path, '*')):
            if aFile in worth_storing:
                shutil.copy(aFile, results_folder_path)
        # run iSS
        iSS_folder_path = path.join(iSS_path, 'results')
        if path.exists(iSS_folder_path):
            shutil.rmtree(iSS_folder_path)
        makedirs(iSS_folder_path)
        output_file = 'OSCAR.DAT'
        if path.isfile(path.join(iSS_path, output_file)):
            remove(path.join(iSS_path, output_file))
        for aFile in glob(path.join(iS_results_path, '*')):
            if aFile in worth_storing:
                shutil.copy(aFile, iSS_folder_path)
        print "running : iSS.e" 
        sys.stdout.flush()
        p = subprocess.Popen('ulimit -n 1000; ./iSS.e', shell=True,
                             stdout=run_record, stderr=err_record, cwd=iSS_path)
        p.wait()
        # backup OSCAR file
        results_folder_name = results_folder_path.split('/')[-1]
        backupiSSOSCAR(iSS_path, results_folder_name)

    # run hydro for v3
    flow_order = 3
    result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_v%d'
                     % (model, ecm, vis, chosen_centrality, tdec, tau0, VisBulkNorm, 
                        eos_name, flow_order))
    if path.isfile(path.join(project_directory, 
                            '%s_%d'%(model, pre_eq),
                            'iSS_backup', '%s.zip'%result_folder)):
        print "hydro search results have already been found, skipping v3..."
    else:    
        print "Start to search v3:"
        hydro_initial_path = path.join(hydro_path, 'Initial')
        cleanUpFolder(hydro_initial_path)
        # run pre-equilibrium
        if(pre_eq == True):
            run_pre_eq(initial_path, chosen_centrality, run_record, err_record, tau0, flow_order)
            for aFile in glob(path.join(rootDir, 'fs/data/result/event_1/%g'%tau0, '*')):
                shutil.move(aFile, hydro_initial_path)
        else:
            shutil.copyfile('%s/sdAvg_order_%d_C%s.dat' % (initial_path, flow_order, chosen_centrality),
                        path.join(hydro_path, 'Initial', 'InitialSd.dat'))

        # hydro
        hydro_results_path = path.join(hydro_path, 'results')
        if path.exists(path.join(hydro_results_path)):
            shutil.rmtree(hydro_results_path)
        makedirs(hydro_results_path)
        cmd = './VISHNew.e'
        if(pre_eq == True):
            args = (' IINIT=2 IEOS=7 iEin=0 iLS=200'
                    + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d visbulknorm=%.6f'
                    % (tau0, edec, vis, norm_factor, pre_eq, VisBulkNorm))
        else:
            args = (' IINIT=2 IEOS=7 iEin=1 iLS=200'
                    + ' T0=%6.4f Edec=%7.5f vis=%6.4f factor=%11.9f initialUread=%d visbulknorm=%.6f'
                    % (tau0, edec, vis, norm_factor, pre_eq, VisBulkNorm))

        print "%s : %s" % (chosen_centrality, cmd + args)
        sys.stdout.flush()
        run_record.write(cmd + args)
        p = subprocess.Popen(cmd + args, shell=True, stdout=run_record,
                             stderr=err_record, cwd=hydro_path)
        p.wait()

        # prepare data backup folder for v3
        results_folder_path = path.join(rootDir, 'RESULTS', result_folder)
        if path.exists(results_folder_path):
            shutil.rmtree(results_folder_path)
        makedirs(results_folder_path)
        # save v3 data
        worth_storing = []
        for aGlob in ['surface.dat', 'dec*.dat', 'ecc*.dat', 'VISH2p1_tec.dat']:
            worth_storing.extend(glob(path.join(hydro_results_path, aGlob)))
        for aFile in glob(path.join(hydro_results_path, '*')):
            if aFile in worth_storing:
                shutil.copy(aFile, results_folder_path)
        # run iSS
        iSS_folder_path = path.join(iSS_path, 'results')
        if path.exists(iSS_folder_path):
            shutil.rmtree(iSS_folder_path)
        makedirs(iSS_folder_path)
        output_file = 'OSCAR.DAT'
        if path.isfile(path.join(iSS_path, output_file)):
            remove(path.join(iSS_path, output_file))
        for aFile in glob(path.join(hydro_results_path, '*')):
            if aFile in worth_storing:
                shutil.copy(aFile, iSS_folder_path)
        print "running : iSS.e" 
        sys.stdout.flush()
        p = subprocess.Popen('ulimit -n 1000; ./iSS.e', shell=True,
                             stdout=run_record, stderr=err_record, cwd=iSS_path)
        p.wait()
        # backup OSCAR file
        results_folder_name = results_folder_path.split('/')[-1]
        backupiSSOSCAR(iSS_path, results_folder_name) 
        # clean up 
        run_record.close()
        err_record.close()
        shutil.copy(path.join(rootDir, run_record_file_name), 
            path.join(results_folder_path))
        shutil.copy(path.join(rootDir, err_record_file_name), 
            path.join(results_folder_path))


def run_hybrid_search_precalculated(model, ecm, vis, tdec, 
               tau0, VisBulkNorm, eos_name, chosen_centrality, pre_eq,
               parallel_mode):
    """
    Shell for parameter search mode simulation after VISHNew running. 
    Divert VISHNew input to iSS and start urqmd thereafter.
    This part runs in parallel mode.
    """
    run_record_file_name = 'run_record_usePrecalculated_search.dat'
    err_record_file_name = 'err_record_usePrecalculated_search.dat'
    run_record = open(path.join(rootDir, run_record_file_name), 'w+')
    err_record = open(path.join(rootDir, err_record_file_name), 'w+')

    # run after burner for v2
    for flow_order in [2,3]:
        print "Start to search v%d:"%flow_order
        result_folder = ('%s%.0fVis%gC%sTdec%gTau%gVisBulkNorm%g_%s_v%d'
                         % (model, ecm, vis, chosen_centrality, tdec, tau0, VisBulkNorm,
                            eos_name, flow_order))
        results_folder_path = path.join(rootDir, 'RESULTS', result_folder)

        # check if parameter search result already exists
        backup_file = path.join(project_directory, '%s_%d'%(model, pre_eq),
            result_folder+'.zip')
        if path.isfile(backup_file):
            print 'skip current run because result already exists!\n'
            sys.stdout.flush()
        else:
            if not path.exists(results_folder_path):
                print "run_hybrid_search_precalculated: no folder %s!"%result_folder
                sys.exit(-1)
            # check if pre-calculated particle list data exists
            particleList_folder_path = path.join(project_directory, 
                '%s_%d'%(model, pre_eq), 'particleList_backup')
            particleList_file_path = path.join(particleList_folder_path,
                '%s.zip'%result_folder)
            if path.isfile(particleList_file_path):
                print "skipping UrQMD because particle_list has been found!"
                sys.stdout.flush()
                # copy data from project folder to results folder
                shutil.copy(particleList_file_path, results_folder_path)
                unzip_cmd = 'unzip -qo %s.zip '%result_folder
                subprocess.call(unzip_cmd, shell=True, cwd = results_folder_path)
                remove(path.join(results_folder_path, '%s.zip'%result_folder))
                # collect db and backup v2 search result
                collectDB_start = time.time()
                collectObservables(result_folder, parallel_mode)
                collectDB_end = time.time()
                print "Collect database running time: %f seconds!"%(collectDB_end-collectDB_start)
                savedFlag = moveDB(result_folder, model, pre_eq)
                if savedFlag==True:
                    remove(particleList_file_path)
            else:
                print "No particle_list file found: %s"%particleList_file_path
                sys.stdout.flush()
                oscar_backup_path, oscar_tmp_path = run_afterBurner(results_folder_path, 
                    chosen_centrality, run_record, err_record,
                    results_folder_path, parallel_mode)
                remove(oscar_backup_path)
                savedFlag = move_particleList(results_folder_path, particleList_folder_path)
                # move iSS file back if particle_list is not saved
                if savedFlag == False:
                    print "particle_list backup failed! move oscar file back!"
                    shutil.copy(oscar_tmp_path, oscar_backup_path.split("/")[0:-1])
                    remove(oscar_tmp_path)

    run_record.close()
    err_record.close()
    shutil.copy(path.join(rootDir, run_record_file_name), 
        path.join(results_folder_path))
    shutil.copy(path.join(rootDir, err_record_file_name), 
        path.join(results_folder_path))


def set_eos(eos_name, tdec):
    """
    This function replace the EOS for the whole simulation
    :param eos_name: the name of EOS
    :param tdec: the decoupling temperature (GeV)
    :return edec: the decoupling energy density (GeV/fm^3) corresponds to tdec
    """
    superMC_eos_path = path.join(rootDir, 'superMC/s95p-PCE')
    fs_eos_path = path.join(rootDir, 'fs/EOS')
    hydro_eos_path = path.join(rootDir, 'VISHNew/EOS/EOS_tables')
    iS_eos_path = path.join(rootDir, 'iS/EOS')
    iSS_eos_path = path.join(rootDir, 'iSS/EOS')
    if eos_name == 's95p-v0-PCE165':
        eos_files_path = path.join(rootDir, 'EOS/EOS_s95p/s95p_convertedtables/s95p-PCE165-v0')
    elif eos_name == 's95p-v1-PCE150':
        eos_files_path = path.join(rootDir, 'EOS/EOS_s95p/s95p_convertedtables/s95p-PCE-v1')
    elif eos_name == 's95p-v1':
        eos_files_path = path.join(rootDir, 'EOS/EOS_s95p/s95p_convertedtables/s95p-v1')
    elif eos_name == 'SM-EOS-Q':
        eos_files_path = path.join(rootDir, 'EOS/SMEOSQ')
    elif eos_name == 'blended':
        eos_files_path = path.join(rootDir, 'EOS/EOS_s95p/s95p_convertedtables/blended')
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



def run_afterBurner(input_folder, cen_string, run_record, err_record, results_folder_path,
                    parallel_mode):
    """
    run iSS + osc2u + urqmd once iS or hydro is finished.
    """
    iSS_path = path.join(rootDir, 'iSS')
    iSS_folder_path = path.join(iSS_path, 'results')
    if path.exists(iSS_folder_path):
        shutil.rmtree(iSS_folder_path)
    makedirs(iSS_folder_path)

    # skip iSS if OSCAR data exists
    iss_backup_path = path.join(project_directory,
        '%s_%d'%(model, pre_eq), 'iSS_backup')
    results_folder_name = results_folder_path.split('/')[-1]
    iss_backup_file_path = path.join(iss_backup_path, "%s.zip"%results_folder_name)
    # location list
    iss_files_location_list = [None, None]
    if path.isfile(iss_backup_file_path):
        iss_files_location_list[0] = iss_backup_file_path
        print "OSCAR file exists! Skipping iSS ... "
        shutil.copy(iss_backup_file_path, iSS_path)
        unzip_cmd = 'unzip -qo %s.zip -d ./'%results_folder_name
        p = subprocess.Popen(unzip_cmd, shell=True,
                             stdout=run_record, stderr=err_record, 
                             cwd=iSS_path)
        p.wait()
        iss_files_location_list[1] = path.join(iSS_path, '%s.zip'%results_folder_name)
    else:
        # copy necessary files to iSS folder
        print "OSCAR file does not exist in: %s"%iss_backup_file_path
        worth_moving = []
        for aGlob in ['surface.dat', 'dec*.dat']:
            worth_moving.extend(glob(path.join(input_folder, aGlob)))
        for aFile in glob(path.join(input_folder, '*')):
            if aFile in worth_moving:
                shutil.copy(aFile, iSS_folder_path) 
        output_file = 'OSCAR.DAT'
        if path.isfile(path.join(iSS_path, output_file)):
            remove(path.join(iSS_path, output_file))
        print "%s : %s" % (cen_string, 'iSS.e')
        sys.stdout.flush()
        iss_start = time.time()
        p = subprocess.Popen('ulimit -n 1000; ./iSS.e', shell=True,
                             stdout=run_record, stderr=err_record, cwd=iSS_path)
        p.wait()
        iss_end = time.time()
        print "iSS run time %f seconds"%(iss_end - iss_start)

        worth_storing = []
        for aGlob in ['*vn*.dat']:
            worth_storing.extend(glob(path.join(iSS_folder_path, aGlob)))
        for aFile in glob(path.join(iSS_folder_path, '*')):
            if aFile in worth_storing:
                shutil.copy(aFile, results_folder_path)
        shutil.rmtree(iSS_folder_path)  # clean up

        # # backup OSCAR file
        # results_folder_name = results_folder_path.split('/')[-1]
        # backupiSSOSCAR(iSS_path, results_folder_name)
    sys.stdout.flush()

    if parallel_mode != 0:
        # run subsequent programs in parallel
        result_files = split_iSS_events(number_of_split = parallel_mode,
                                        output_folder = results_folder_path)
    else:
        #osc2u
        o2u_path = path.join(rootDir,'osc2u')
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
        UrQMD_path = path.abspath(rootDir,'urqmd')
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
        urqmd_start = time.time()
        p = subprocess.Popen('bash runqmd.sh', shell=True, stdout=run_record,
                             stderr=err_record, cwd=UrQMD_path)
        p.wait()
        urqmd_end = time.time()
        print "UrQMD running time: %f seconds!"%(urqmd_end-urqmd_start)

        worth_storing = []
        for aGlob in ['particle_list.dat']:
            worth_storing.extend(glob(path.join(UrQMD_path, aGlob)))
        for aFile in glob(path.join(UrQMD_path, '*')):
            if aFile in worth_storing:
                shutil.copy(aFile, results_folder_path)
        remove(path.join(UrQMD_path, input_file))  # clean up
        remove(path.join(UrQMD_path, output_file))  # clean up
    return iss_files_location_list

def move_particleList(source_folder_path, backup_folder_path):
    """
    zip and move urqmd generated particle_list files to the project folder
    """
    particle_list_pattern = 'particle_list*.dat'
    results_folder_name = source_folder_path.split('/')[-1]
    # zip the particle_list file
    zip_cmd = 'zip -q %s.zip %s'%(results_folder_name, particle_list_pattern)
    subprocess.call(zip_cmd, shell=True, cwd = source_folder_path)
    # move the file
    if path.exists(path.join(backup_folder_path, '%s.zip'%results_folder_name)):
        remove(path.join(backup_folder_path, '%s.zip'%results_folder_name))
    shutil.move(path.join(source_folder_path, '%s.zip'%results_folder_name),
                backup_folder_path)
    if path.isfile(path.join(backup_folder_path, '%s.zip'%results_folder_name)):
        savedFlag = True
        print "urqmd particle_list file saved!"
    else:
        print "urqmd particle_list file move failed!"
        savedFlag = False
    return savedFlag

def run_simulations(mode, model, ecm, dN_deta, vis, tdec, tau0, VisBulkNorm, eos_name,
                    cf_flag, fit_flag, chosen_centrality, collsys, pre_eq,
                    parallel_mode):
    """
    shell function to run simulations
    :param mode: simulation mode: hydro or hybrid
    :param model: initial condition model
    :param ecm: collision energy
    :param dN_deta: final charged multiplicity
    :param vis: the specific shear viscosity
    :param tdec: the decoupling temperature
    :param tau0: the starting time of hydrodynamic simulation
    :param VisBulkNorm: normalization factor for parametrized zeta/s(T)
    :param eos_name: the name of EOS
    :param cf_flag: switch for Cooper-Frye freeze-out
    :param pre_eq : switch for pre-equilibrium evolution
    :return: none
    """
    print('%s mode: %s sqrt{s} = %s A GeV' % (mode, model, ecm))
    print('eta/s = %g, Tdec = %g GeV, tau0 = %g fm/c' % (vis, tdec, tau0))
    print('zeta/s norm factor = %g'%VisBulkNorm)
    print('Pre-equilibrium: %s'%pre_eq)
    print('EOS : %s' % eos_name)

    # initial setup
    result_folder_path = path.join(rootDir, 'RESULTS')
    if not path.exists(result_folder_path):
        makedirs(result_folder_path)
    init_conditions_path = path.join(result_folder_path, 'initial_conditions')
    if not path.exists(init_conditions_path):
        makedirs(init_conditions_path)

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
                                           collsys, pre_eq=pre_eq)
        else:
            sys.exit(0)
    else:
        if modelsys == "MCKLN":
            preEq_string = "_preEq" if pre_eq == True else "_nopreEq"
            initial_condition_name = '%s%.0f_sigmaNN_gauss_d0.9%s' % (modelsys, ecm, preEq_string)
        else:
            initial_condition_name = '%s%.0f_sigmaNN_gauss_d0.9' % (modelsys, ecm)            
        if path.isfile(path.join(rootDir, 'initial_conditions/%s.zip' % initial_condition_name)):
            p = subprocess.Popen('unzip -qo %s.zip' % initial_condition_name, 
                                 shell=True, stdout=subprocess.PIPE, 
                                 cwd=path.join(rootDir, 'initial_conditions'))
            p.wait()
            # move initial profiles
            move_cmd = "mv -f ./initial_conditions/%s/* %s/initial_conditions/"%(
                initial_condition_name, result_folder_path)
            p = subprocess.Popen(move_cmd, shell=True, stdout=subprocess.PIPE,
                cwd=rootDir)
            p.wait()
        else:
            print("initial density profiles for %s%% centrality for %s %s " 
                  "at sqrt{s} = %g A GeV is not found!" 
                  % (chosen_centrality, model, collsys, ecm))
            generate_flag = 'y'#raw_input("Do you want to generate one right now?")
            if generate_flag.lower() in ['yes', 'y']:
                generate_avg_initial_condition(model, ecm, chosen_centrality, 
                                               collsys, pre_eq=pre_eq)
            else:
                sys.exit(0)

    # start to run simulations
    if fit_flag:
        print "fitting the overall normalization factor ..."
        norm_factor_guess = linearFitSfactor(tau0, vis, tdec, model, pre_eq) # guess the scaling factor from linear regression
        if norm_factor_guess <=0:
            print "predicted sfactor smaller than 0: "+"sfactor=%g"%norm_factor_guess
            norm_factor_guess = norm_factor_default
        if mode == "hybrid_search" or mode == "hydro_search":
            norm_factor = fit_hydro(dN_deta, vis, edec, tau0, VisBulkNorm, pre_eq, norm_factor_guess, chosen_centrality) # directly fit hydro at 10-20% centrality
        else:
            norm_factor = fit_hydro(dN_deta, vis, edec, tau0, VisBulkNorm, pre_eq, norm_factor_guess)
    else:
        norm_factor = norm_factor_default
    if mode == 'hydro':
        print "running pure hydro simulations for centrality bin(s): %s ..."%chosen_centrality
        run_purehydro(model, ecm, norm_factor, vis, tdec, edec, tau0, VisBulkNorm,
                      eos_name, cf_flag, chosen_centrality, pre_eq)
    elif mode == 'hybrid':
        print "running hybrid simulations for centrality bin(s): %s..."%chosen_centrality
        run_hybrid(model, ecm, norm_factor, vis, tdec, edec, tau0, VisBulkNorm,
                   eos_name, chosen_centrality, pre_eq, parallel_mode)
    elif mode =="hybrid_search":
        print "running hybrid search simulations for centrality bin(s): %s..."%chosen_centrality
        run_hybrid_search(model, ecm, norm_factor, vis, tdec, edec, tau0, VisBulkNorm,
                   eos_name, chosen_centrality, pre_eq, parallel_mode)
    elif mode =="hydro_search":
        print "running pure hydro search simulations for centrality bin(s): %s..."%chosen_centrality
        run_hydro_search(model, ecm, norm_factor, vis, tdec, edec, tau0, VisBulkNorm,
                   eos_name, chosen_centrality, pre_eq)        
    else:
        print sys.argv[0], ': invalid running mode', mode
        sys.exit(1)


def print_help_message():
    print "Usage : "
    print(color.bold
          + "./runHydro.py -ecm ecm "
          + "[-mode mode -model model -vis vis -Tdec Tdec -tau0 tau0 -VisBulkNorm VisBulkNorm"
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
    print(color.bold + "-VisBulkNorm" + color.end
          + "  the normalization factor for zeta/s(T) \n"
          + color.bold + "       VisBulkNorm = 1.0 [default]" + color.end)
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
    print(color.bold + "-parallel_mode" + color.end
          + "   switch to split iSS events and run osc2u+urqmd in parallel: "
          + color.purple + color.bold + "0 [default]" + color.end
          + color.purple + ', e.g. integer N (>1) to split iSS events to N pieces' + color.end)
    print(color.bold + "-h | -help" + color.end + "    This message")


if __name__ == "__main__":
    # set default values
    vis = 0.08
    tdec = 0.12  # GeV
    tau0 = 0.6  # fm/c
    VisBulkNorm = 1.0 
    mode = 'hydro'
    model = 'MCGlb'
    eos_name = 's95p-v0-PCE165'
    chosen_centrality = 'All'
    collsys = 'Pb+Pb'
    cf_flag = True
    fit_flag = True
    pre_eq = False
    parallel_mode = 0

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
        elif option == '-VisBulkNorm':
            VisBulkNorm = float(sys.argv[1])
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
        elif option == '-parallel_mode':
            parallel_mode  = int(sys.argv[1])
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
    norm_factor_default = 56.76 if model == 'MCGlb' else 9.92

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
        if mode=='hybrid_search' or mode =='hydro_search':
            if ecm_string=='2760.0':
                dN_deta = dn_deta_ECM2760_dict[chosen_centrality]
            else:
                print 'Mode hybrid_search/hydro_search is currently not supported for %s GeV'%ecm_string
                sys.exit(0)
    else:
        print sys.argv[0], ': invalid collision energy', ecm
        sys.exit(1)

    if mode in ['hydro', 'hybrid', 'hybrid_search', 'hydro_search']:
        run_simulations(mode, model, ecm, dN_deta, vis, tdec, tau0, VisBulkNorm, eos_name,
                        cf_flag, fit_flag, chosen_centrality, collsys, pre_eq,
                        parallel_mode)
    elif mode == 'hybrid_usePrecalculated':
        run_hybrid_search_precalculated(model, ecm, vis, tdec, 
               tau0, VisBulkNorm, eos_name, chosen_centrality, pre_eq,
               parallel_mode)
    else:
        print sys.argv[0], ': invalid running mode', mode
        print_help_message()
        sys.exit(1)
