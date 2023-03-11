#!/usr/bin/env python3
def select_funct(str_ex,str_ec) :
    """ This function selects exchange and correlation functionals into 
    FLOSIC/NRLMOL format """
    funct = str_ex+str_ec
    str1 = str_ex.lower()
    str2 = str_ec.lower() 
    if str1 == 'scan' :
        funct = 'MGGA-SCAN'

    elif str1 == 'rscan' :
        funct = 'MGGA-RSCAN'

    elif str1 == 'r2scan' :
        funct = 'MGGA-R2SCAN'

    elif str1 == 'lda' :
        funct = 'LDA-PW91'

    elif str1 == 'pbeibp' :
        funct = 'MGGA-PBEIBP'

    elif str1 ==  'pbe' :
        funct = 'GGA-PBE'

    funct = funct + '*'

    if str2 == 'scan' :
        funct = funct + 'MGGA-SCAN'

    elif str2 == 'rscan' :
        funct = funct + 'MGGA-RSCAN'

    elif str2 == 'r2scan' :
        funct = funct + 'MGGA-R2SCAN'

    elif str2 == 'lda' :
        funct = funct + 'LDA-PW91'

    elif str2 == 'pbeibp' :
        funct = funct + 'MGGA-PBEIBP'

    elif str2 ==  'pbe' :
        funct = funct + 'GGA-PBE'

    return funct 


def get_binary_name(code) :
    """ This function provides the binary name of the version code  """ 
    bn=''
    code1 = code.lower()
    if code1 == 'flosic1' :
        bn = 'nrlmol_exe'

    if code1 == 'flosic2' :
        bn = 'flosic_exe'

    if code1 == 'blockmesh' :
        bn = 'dftsw_exe'

    return bn

def selec_set(test_set) :
    """ This function provides the dataset to run according to argument sent"""
    set1=test_set.lower()
    arr = []
    if set1 == 'atoms' :
        for i in range(1,37) :
            arr.append(i)

    elif set1 == 'ae6' :
        for i in range(128,134) :
            arr.append(i)

    elif set1 == 'anions' :
        arr = [1, 3, 5, 6, 8, 9, 11, 13, 14, 15, 16, 17, 19, 22, 29, 31, 32, 33, 34, 35]

    elif set1 == 'cations' :
        arr = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, \
                 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
 
    elif set1 == 'molec' :
        for i in range(101,146) :
            arr.append(i)
 
    elif set1 == 'bh6' :
        arr = ['CH3','CH4', 'H2','H2O','H2S','HS','OH','saddle1','saddle2','saddle3'] 


    elif set1 == 'noble' :
        arr = ['He.in','Ne.in','Ar.in','Kr.in','Xe.in']

    
    elif set1 == 'wcpt18' :
        arr = [ 'h2o','reac1','reac2','reac3','reac4','reac5','reac6','reac7','reac8','reac9', \
                'ts1','ts1h2o','ts2','ts2h2o','ts3','ts3h2o','ts4','ts4h2o','ts5','ts5h2o','ts6',\
                'ts6h2o','ts7','ts7h2o','ts8','ts8h2o','ts9','ts9h2o' ] 
 
    elif set1 == 'sie' :
        for i in range(1,37) :
            arr.append(i)
        arr.append = ['H2O','H2Op','He','Hep','NH3','NH3p']

    elif set1 == 'gw100' :
        for i in range(1,101) :
            arr.append(i)

    return arr

def execute_bin( server, path1, subset, niter, nproc) :
    """ This function will run the binary in path1 according to server
        server : linux, perlmutter or lonestar"""
    n = len(niter)

    if server == 'linux' : 
        for irun in niter : 
            exect = 'mpirun -np '
            exect = exect + str(nproc) +'  '+ path1 + ' > print.'+str(irun)
            os.system(exect) 

    elif server == 'perlmutter' :
            nproc = 128
            jobf='job_perlmutter.sl'
            f = open(jobf,'w') 
            f.write("#!/bin/bash \n")
            f.write("#SBATCH -A m3714 \n")
            f.write("#SBATCH -C cpu \n")
            f.write("#SBATCH -N 1 \n")
            f.write("#SBATCH -q regular \n")
            f.write("#SBATCH -J " + subset + "\n" )
            f.write("#SBATCH --mail-user=ssromerogon@miners.utep.edu \n")
            f.write("#SBATCH -t 12:00:00 \n\n")
            f.write("#OpenMP settings: \n")
            f.write("export OMP_NUM_THREADS=1 \n")
            f.write("export OMP_PLACES=threads  \n")
            f.write("export OMP_PROC_BIND=spread \n")
            f.write('export SLURM_CPU_BIND="cores" \n\n')
            f.write("#run the application:    \n")
            f.write("for i in {1.." + str(n) + "}; do \n")
            f.write('   if [ -e "EXIT" ]; then \n')
            f.write("      rm -f EXIT; exit    \n")
            f.write("   fi                     \n")
            f.write('   echo "1" > igroup      \n')
            f.write('   echo "'+str(nproc)+'" >> igroup   \n')
            exect = "   srun -n " + str(nproc) + " -c 2 --cpu_bind=cores /global/homes/s" + path1 + " > print.$i \n"
            f.write(exect)
            f.write('   if [ -e "NRLMOL_ERROR" ]; then \n')
            f.write("      rm -f NRLMOL_ERROR; exit    \n")
            f.write("   fi                    \n")
            f.write('   if [ -e "DFTSW_ERROR" ]; then \n')
            f.write("      rm -f DFTSW_ERROR; exit   \n")
            f.write("   fi                    \n")
            f.write("done \n") 
            f.close()
            exect = 'sbatch ' + jobf
            os.system(exect)

    return True

def create_input( blend, guess_fod, fod_loop, fix1s, max_scf, basis, calct, scf_tol, tmesh, code) :
    """ This function creates a default NRLMOL/DFTSW.INP file """
    code1 = code.lower()
    if (code1 == 'flosic1') or (code1 == 'flosic2') :
            f = open('NRLMOL_INPUT.DAT','w')
            f.write("&input_data\n")
            f.write("ATCALCV       = '-1' ! list of functionals for @Calculations (-1 to bypass)\n")
            f.write("ATOMSPHV      = 'N'\n")
            f.write("BASISV        = "+basis+" ! Specify basis for calculation(basis.txt)\n")
            f.write("BLENDERV      = "+blend+" ! FOD blender mode (SIC only)\n")
            f.write("CALCTYPEV     = "+calct+"\n")
            f.write("DFTD3V        = 'N' ! Set to Y to do include Grimmes DFT-D3 dispersion\n")
            f.write("DIAG1V        =  1  ! diagonalization to use on regular arrays (diagge.f90)\n")
            f.write("DIAG2V        =  1  ! diagonalization to use on packed arrays (diag_dspgv.f90)\n")
            f.write("DIAG3V        =  0  ! diagonalization to use on parallel (sdiagge_n.f90)\n")
            f.write("DMATV         = 'N' ! Create/use/mix density matrix\n")
            f.write("DOSOCCUV      = 'N' ! Controls wether to calculate density of states\n")
            f.write("EXCITEDV      = 'N' ! Determines if this is an excited state calculation (DFA only)\n")
            f.write("FIXMV         = 'Y' ! Fix spin moment\n")
            f.write("FORMFAKV      = 'N' ! this controls if FORMFAK is executed\n")
            f.write("FOD_LOOPV     = "+fod_loop+" ! Inernal FOD loop for frozen density optimization\n")
            f.write("FOD_OPT1V     = 'LBFGS' ! FOD_OPT: algorithm (LBFGS/CG)\n")
            f.write("FOD_OPT2V     = 'N' ! FOD_OPT: scaling of r and F\n")
            f.write("FOD_OPT3V     = "+str(fix1s)+"  ! FOD_OPT (0)no constraint (1)fix1s (2)fullForce (3)shellOPT (4)freeFOD\n")
            f.write("FRAGMENTV     = 'N' ! Process CLUSTER in fragments\n")
            f.write("FROZENDENV    = 'N' ! Frozen density mode (SIC only)\n")
            f.write("JNTDOSV       = 'N' ! This calculates joint density of states (DFA only)\n")
            f.write("MATDIPOLEV    = 'N'\n")
            f.write("MAXSCFV       = "+str(max_scf)+" ! Maximum SCF iterations\n")
            f.write("MIXINGV       = 'P' ! (H)amiltonian (P)otential (D)ensity matrix mixing\n")
            f.write("MOLDENV       = 'N' ! Use molden and wfx driver\n")
            f.write("MPI_IOV       = 'N' ! Use MPI IO for scalapack distribution\n")
            f.write("NBOV          = 'N' ! Use NBO driver\n")
            f.write("NONSCFV       = 'N' ! Set to Y to do a non SCF calculation\n")
            f.write("NONSCFFORCESV = 'N' ! Set to Y to calculate forces in a non SCF calculation\n")
            f.write("NWFOUTV       = 10    ! Write WFOUT file for every N-th iteration\n")
            f.write("POPULATIONV   = 'N' ! Population analysis\n")
            f.write("RHOGRIDV      = 'N' ! Set to Y to execute RHOGRID\n")
            f.write("SCFTOLV       = "+scf_tol+" ! SCF tolerance\n")
            f.write("SPNORBV       = 'N' ! Run SPNORB (DFA only)\n")
            f.write("SPNPOLV       = 'N' ! Run spin polarized calculation from CLUSTER\n")
            f.write("SOLVENTV      = 'N' ! Set to Y to include solvent effect (SOLVENTS)\n")
            f.write("SCANMESHV     = "+str(tmesh)+"   ! Mesh recommended for (0)LDA/PBE, (1)SCAN (2)rSCAN\n")
            f.write("SYMMETRYV     = 'N' ! Set to Y to detect symmetry\n")
            f.write("UNIHAMV       = 'N' ! Set to Y to use unified Hamiltonian formalism in SCF-SIC (SIC only)\n")
            f.write("WFGRIDV       = 'N' ! set to Y to write orbitals in cube format (DFA only)\n")
            f.write("WFFRMV        = 'N' ! set to Y to write Fermi orbitals in cube format (SIC only)\n")
            f.write("GUESS_FODV    = "+guess_fod+" ! set to Y to get a guess FOD set\n")
            f.write("ADSICSTARTV   = 'N' ! Set to Y for ADSIC starting guess (SIC only)\n")
            f.write("&end\n")
            f.close()

    #if code1 == 'blockmesh'      

    return True


def format_str(lines) :
    '''This function removes spaces and makes it tuple '''
    new_format = []
    for line in lines :
        #ln = [ x.strip() for x in line.split(',')]
        ln = re.split(r'\s+', line)
        new_format.append(ln)

    return new_format


def check_fod_conv(parse) :  
    """ This function will check if FODs have converged by checking fande.out with enough lenght """
    #Check if fande.out exists, if not just convergence false  
    file1 = open('fande.out', 'r')
    Lines = file1.readlines()
    #Parsing Evalues? 
    parse = False
    
    results = []
    
    count = 0
    # Strips the newline character, reading fande.out into  results as string tuple 
    for line in Lines:
        count += 1
        results.append((line.strip()))
 
    # Print fande.dat if parse is True  
    if parse : 
        for line in results:
            print(line)
    # Converting results to colum separated tuple 
    res_formatted = format_str(results)
    if parse : 
        for line in res_formatted :
            print(line)
    
    fod_force_cols = len(res_formatted[0][:])

    #print("FOD file contains this number of columns :  ", fod_force_cols )
    # Passing the results to a dataframe
    if fod_force_cols == 3 :
        fod_df = pd.DataFrame(res_formatted, columns = ['Iter','Ener','max_force']) 
        fod_df = fod_df.astype({'Iter':'int','Ener':'float','max_force':'float'})
    else :  
        fod_df = pd.DataFrame(res_formatted, columns = ['Iter','Ener','sqr_force','max_force']) 
        fod_df = fod_df.astype({'Iter':'int','Ener':'float','sqr_force':'float','max_force':'float'})

    # This cleans inconsistency within frozen density ( 3 cols ) and regular fod output (4 cols)
    #fod_df['max_force'] = fod_df['max_force'].fillna(-1000.0)
    
    niter = len(fod_df.loc[:,"Iter"]) 
    #if niter == 1 :
    #    print("First FOD update, not possible to determine...") 
    
    fod_energy = 10000.0
    fod_force = 10000.0
    last_iter = -1
    max_force = -1000.0
    conv = False
    Ediff_criteria = 0.000005
    Fdiff_criteria = 0.0005
    fod_criteria = 0.002
    fod_full_criteria = 0.0005
    str2 = ''
    str2_old = ''
    # Tight criteria convergence (this worked for CO spin split) 
    for i in range( 1, niter) :
        str1 = 'Ener'
        # Compatibility with frozen density of regular optimization
        if fod_df.loc[i,'max_force'] == -1000.0  :
            str2 = 'sqr_force'
            if i <= 1 : 
                str2_old = str2
        else :
            str2 = 'max_force'
            if i <= 1 : 
                str2_old = str2
     
        Ediff = abs( fod_df.loc[i, str1] - fod_df.loc[i-1, str1])
        Fdiff = abs( fod_df.loc[i, str2] - fod_df.loc[i-1, str2])
        max_force = abs( fod_df.loc[ i, str2])
        #print('Vals ans Current differences : ', fod_df.loc[i,str1],fod_df.loc[i,str2],Ediff, Fdiff ) 
        # This three criteria together are quite tight
        if ( ((Ediff <= Ediff_criteria) and \
              (Fdiff <= Fdiff_criteria) and \
              (max_force <= fod_criteria ) ) or \
              (max_force <= fod_full_criteria) ) :
            f=open('FOD_CONV','w')
            f.write("Converged:\t {} \n".format(i))
            f.write("Ener&fforce:\t {} \t {} \n".format(fod_df.loc[i,str1],fod_df.loc[i,str2]) )
            f.write("E&Fdiff:\t {} \t {} \n".format(Ediff, Fdiff) )
            f.close()
            fod_energy = float(fod_df.loc[i,str1])
            fod_force = float(fod_df.loc[i,str2])
            last_iter = i
            conv = True
            print('Energy and max_fod_force', fod_energy, fod_force) 
            break 
    
    if conv :
        print('FODs have converged')
    else :
        print("FODs not converged yet... cheer up" )
    
    fod_conv_res = [fod_energy, fod_force, conv, last_iter]
    return fod_conv_res 

def get_evalues(parse) :
    """ This function will retrieve eigenvalues (HOMO and LUMO) from EVALUES file  """
    file1 = open('EVALUES', 'r')
    Lines = file1.readlines()
    #Parsing Evalues? 
    
    
    str2search = 'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
    results = []
    
    count = 0
    ltest = False
    # Strips the newline character
    for line in Lines:
        count += 1
        if ltest :
            results.append((line.strip())) 
        if str2search in line:
            ltest=True
            #print("Line{}: {}".format(count, line.strip()))
    if parse : 
        for line in results:
            print(line)
    
     
    #Reformat results string to tuples  
    res_formatted = format_str(results)

    orb_info = []
    occu = []
    count=0
    occu=0
    homo=0
    lumo=0 
    for line in results:
        # Store orbital, spin, energy, occupation
        orb_info.append( ( int( res_formatted[count][0] ) , int(res_formatted[count][6]), float(res_formatted[count][8]), float(res_formatted[count][10]) ) ) 
        #print(orb_info[count], 'Type : ', type(orb_info[count]) )  
        occu += orb_info[count][3]
        if orb_info[count][3] == 0.0 :
            homo = orb_info[count-1][2] 
            lumo = orb_info[count][2] 
            break;
        count += 1
    
    print('HOMO&LUMO:', homo,lumo)
    evals=[homo,lumo]

    return evals    

#def get_energy(parse) :
#    """ This function will retrieve total energy, sic energy and dfa energy """
#    file1 = open('SUMMARY', 'r')
#    Lines = file1.readlines()
#    #Parsing Evalues? 
#    
#    
#    results = []
#    
#    count = 0
#    ltest = False
#    # Strips the newline character
#    for line in Lines:
#        results.append((line.strip())) 
#        count += 1
#    if parse : 
#        for line in results:
#            print(line)
#    
#     
#    #Reformat results string to tuples  
#    res_formatted = format_str(results)
#
#    ener_info = []
#    count=0
#    for line in results:
#        # Store orbital, spin, energy, occupation
#        if not res_formatted[count][0] == 'IT' :
#            ener_info.append( ( float( res_formatted[count][2] ) , float(res_formatted[count][5]) ) ) 
#            print(ener_info[count], 'Type : ', type(ener_info[count]) )  
#        
#        #print(res_formatted[count]) 
#        count += 1
#    
#    # Do energy convrgence later ... 
#    # this is true for FLOSIC, not for DFTSW
#    Edfa = ener_info[icount-1][0]
#    Etot = ener_info[icount-1][5]
#    Esic = Etot - Edfa
#    eners=[ Etot, Esic, Edfa]
#    print("Last energy results : ",eners ) 
#    exit()
#    return eners    




import pandas as pd
import numpy as np
import os, sys, math, shutil, re

###################################################
############# Enviroment setup  #####################
server   = 'linux'    #perlmutter, ls6, linux
username = 'sromero'
# Will clean the running path 
clean    = False 
# Write results only and skip computation
res_only = False
nproc    = 4 
niter    = range(1,2) # Run this many times the binary 
code     = 'flosic2'  # flosic1 = nrlmol_exe, flosic2 = flosic_exe  or blockmesh = dftsw_exe
subdirs  = ['base']
binary_loc = '/0.3'   # from /home/username
##################################################################
######################### Running setup #########################
# Basis set for FLOSIC or BLOCKMESH
basis   = "'DEFAULT'" # set it to write
# Exchange and correlation
exch    = 'lda'
corr    = 'lda'
max_scf = 200
calct   = "'SCF-ONLY'" # can also be LBFGS for molecular optimization
scf_tol = "1.0D-6" 
tmesh   = 0 # 0 regular mesh, 1 scan mesh, 2 rscan mesh
# FOD setup
blend   = "'Y'"
fix1s   = 1 # 0 all electrons, 1 fix1s
guess_fod   = "'N'" 
fod_loop    = "'N'" 
##################################################################


# Get info 
narg        = len(sys.argv)
test_set    = sys.argv[1]
# Set must come from first argument 
set_list    = selec_set(test_set)
# Absolute path for binary
bpath       ='/home/' + username + binary_loc
bname       = get_binary_name(code)
# functional name defined here
funct       = select_funct(exch,corr)

# Custom binaries for different subdirs  
# bname = ['flosic_exe', 'flosic_exe_100', 'flosic_70' ]
# subdirs  = ['sic', 'sosic_100', 'sosic_70' ]

# Check if set_list is void and exit
if len(set_list) == 0 :
    print('Benchmark set {} not included, include it by hand or use an included set'.format(test_set))
    exit()

print("=================================================================")
print("Running info ...")
print("Functional  : ", funct)
print("Binary path : ", bpath, " Binary name : ",bname )
print('Set : ', test_set ,set_list)
if clean : 
    print(" WARNING : cleaning mode !!!") 
print("=================================================================")

# main running path 
main_path = os.getcwd()
ener_res = []
fod_res = []
evals_res = []

for jcalc in set_list :

    os.chdir( str(jcalc) )
    curr_path = os.getcwd() 

    # Copy to input directory all initial files of benchmark
    if not os.path.exists( 'input' ) :

        files = os.listdir( curr_path )
        os.makedirs( 'input' )
        # Move files from jcalc curr_path to input subdirectory
        for file in files : 
            target_dir = os.path.join('input/', file)
            shutil.move(file,target_dir) 

    if clean :
        for dirs in subdirs :
            rmf = os.path.join( curr_path, dirs)
            shutil.rmtree(rmf)

    # Go into sub directories
    for dirs in subdirs :

        # Create subdirectory if needed 
        if not os.path.exists( str(dirs) ) : 
            os.makedirs( str(dirs) )
        # Move to subdirectory path
        os.chdir( str(dirs) )
        curr_path2 = os.getcwd()
        
        # Check if input files are there, if not bring those... 
        files = os.listdir(curr_path2)
        fmatch = 'SYMBOL' 
        xmatch = [ifile for ifile in files if ifile == fmatch ]
        f1exist = any(xmatch)
        if not f1exist :  
            fmatch = 'CLUSTER'
            xmatch = [ifile for ifile in files if ifile == fmatch ]
            f1exist = any(xmatch)
        fmatch = 'SUMMARY' 
        xmatch = [ifile for ifile in files if ifile == fmatch ]
        f2exist = any(xmatch)
        first_run = (not f1exist) or (not f2exist)
        # Not enough, if there is the symbol file we would like to rerun
               
        # Copy files from input directory  : NOTE copy only when wanted
        if first_run :
            esc_calc = False  
            files = os.listdir(curr_path+'/input')
            for file in files : 
                shutil.copy2( curr_path + '/input/'+ file, curr_path2)
        

        # Set exchange correlation
        if os.path.exists('SYMBOL') :
            os.system('sed -i "2s|.*|' + funct + '|" SYMBOL ')
        if os.path.exists('CLUSTER') :
            os.system('sed -i "1s|.*|' + funct + '|" CLUSTER ')

        if first_run : 
            blend = "'N'" 
        if (blend == "'Y'") and (not first_run) : 
            # start from sctratch 
            os.system('sed -i "2s|.*| 3  4 |" RUNS ')
 
        # Create NRLMOL_INPUT.DAT
        create_input( blend, guess_fod, fod_loop, fix1s, max_scf, basis, calct, scf_tol, tmesh, code) 
 
        path1 = bpath + '/' + bname
        # Here will run server = machine
        # path1 is binary path + binary name 
        # jcalc is job name if needed for HPC (perlmutter or lonestar)
        # niter is number of times to run binary
        if not first_run : 
            # FOD convergence criteria
            fod_res1  = check_fod_conv(False)
            fod_res.append(fod_res1)
            # obtain eigenvalues (beta version) 
            evals1 = get_evalues(False)
            evals_res.append(evals1)
            # Obtain last energy in Summary
            #ener_res1 = get_energy(False) 
            #ener_res.append(ener_res1)
            # Obtain dipole moments

            # If FODs converged skip calculation 
            esc_calc = fod_res1[2]

        # Only write results and skip computation    
        if res_only :
            esc_calc = True

        if not esc_calc :
            print("Running ", jcalc, "..." )
            # if read true will skip the running/executing part
            execute_bin(server, path1, str(jcalc), niter, nproc)

        # Go back to jcalc path
        print(" ")
        os.chdir(curr_path)

    os.chdir( main_path )

#ener_res_df = pd.DataFrame( ener_res, columns = ['Etot', 'Edfa', 'Esic']) 
#ener_res_df.to_csv('energy.csv', index=False)
fod_res_df = pd.DataFrame( fod_res, columns = ['Ener', 'max_force', 'conv','conv_iter']) 
fod_res_df.to_csv('fodres.csv', index=True, sep='\t')
evals_res_df = pd.DataFrame( evals_res, columns = ['Homo', 'Lumo'] )
evals_res_df.to_csv('evals.csv', index=True, sep='\t')

