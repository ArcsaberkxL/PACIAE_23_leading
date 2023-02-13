#!/bin/bash
#SBATCH -J PACIAE_SLURM       # Task running name. Arbitrary.
#SBATCH -N 1                  # Number of nodes. Modify as needed.
#SBATCH --ntasks-per-node=1   # Number of cores used per node. Modify as needed.
#SBATCH -p middle             # Name of sequence. Modify as needed.

# Above statements are for normal LINUX system in PC and the SLURM scheduling 
#   system on the computing cluster or super-computer.
# Do not delete them but modify them as needed.

# For LSF scheduling system on the computing clusters or super-computers.
APP_NAME=Gsx_normal  # Name of sequence. Modify as needed.
NP=18                # Total number of cores used. Modify as needed.
NP_PER_NODE=18       # Number of cores used per node. Modify as needed.
RUN="RAW"            # Additional option. Not required to be modified usually.

################################################################################
################################################################################
###                                                                          ###
###   PPPPPPP       AAAAA       CCCCCCC    IIIIIII     AAAAA     EEEEEEEEE   ###
###   P      p     A     A     C       C      I       A     A    E           ###
###   P       p   A       A   C               I      A       A   E           ###
###   P      p    A       A   C               I      A       A   E           ###
###   PPPPPPP     AAAAAAAAA   C               I      AAAAAAAAA   EEEEEEEEE   ###
###   P           A       A   C               I      A       A   E           ###
###   P           A       A   C               I      A       A   E           ###
###   P           A       A    C       C      I      A       A   E           ###
###   P           A       A     CCCCCCC    IIIIIII   A       A   EEEEEEEEE   ###
###                                                                          ###
################################################################################
################################################################################
#                                                                              #
#                                                                              #
################################################################################
#                                                                              #
# This is a universal toy SHELL-script for PACIAE running on normal LINUX      #
#   system in the personal-computer and task submitting on SLURM and LSF       #
#   scheduling systems in the computing cluster and the sunper-computer.       #
#                                                                              #
# Mainly for NN, NA and AA collisions.                                         #
#                                                                              #
# This shell script will generate usu.dat and Makefile automatically. Then     #
#   run / submit tasks automatically, too.                                     #
# A series of folders will be generated to store different type of files,      #
#   including src (source folder, stores the source code files), bin (binary   #
#   folder, stores binary executable files), log (log records folder, stores   #
#   running informations) , and etc (editable text configuration folder,       #
#   stores configuration and other files). The PACIAE running folder with the  #
#   name of collisions system is also generated for easier management.         #
#                                                                              #
# It makes the (pseudo-)SPMD (single program, multiple data) possible to       #
#   achieve the (pseudo-)parallel running.                                     #
# Total N events are split up to multiple separate (N/number_of_CPU) events    #
#   and run simultaneously on multiple processors with different input random  #
#   seed based on the real-time clock of the machine.                          #
# This random number generator seed is set in PACIAE internal code (main.f).   #
#                                                                              #
# One needs to note that it's just a pseudo-parallism. In other words, this    #
#   script just performs PACIAE running one-by-one automatically instead of    #
#   manual runs. One also needs another program to aggregate and average all   #
#   the results (in rms.out) after the whole PACIAE running finished.          #
#                                                                              #
# This script has been tested on UBUNTU 20.04 in PC, on SLURM in CCNU          #
#   FARM-cluster (middle partition), and on LSF in National SuperComputing     #
#   Center in Shenzhen (NSCCSZ, Gsx_normal partition).                         #
#                                                                              #
################################################################################
#                                                                              #
# How to use:                                                                  #
#                                                                              #
#   1. First and formost, give this script file "executable permission" by     #
#        typing command "chmod +x PACIAE.sh". In addition, one needs "make"    #
#        tool. Install "make" using the command, for example on UBUNTU,        #
#        "sudo apt install make".                                              #
#                                                                              #
#   2. Modify the variabls needed.                                             #
#        2.1 Modify the names of source files for Makefile. They shoud be the  #
#              same as the real ones (except for name_x, it's user-defined).   #
#        2.2 Modify the variables: the number of program-running (n_run,       #
#              preferably not larger than the total number of CPU cores of the #
#              computer), the number of events per program-running (n_eve),    #
#              and other main setups of incident particles and collisions etc. #
#              "(D=)" means default value.                                     #
#            Here, the variables are aliases for the ones in usu.dat.          #
#            More detailed setups can be modified in the closely following     #
#              usu.dat directly. Search keywords "usu.dat" to find them.       #
#        2.3 On SLURM/LSF scheduling systems in the computing cluster and the  #
#              sunper-computer. One needs modify the additional statements of  #
#              the settings at the beginning of file. (SBATCH, APP_NAME...)    #
#                                                                              #
#   3. Execution.                                                              #
#        3.1 On normal LINUX, to run this script by typing command             #
#              "./PACIAE.sh".                                                  #
#            The more recommended command is                                   #
#              "time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log",           #
#              which stores the screen information to a log file.              #
#            Also, one can execute them in the backgroud using                 #
#              "nohup ./PACIAE.sh &" or                                        #
#              "nohup time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log &"    #
#       3.2 On SLURM: "sbatch PACIAE.sh". On LSF: "bsub PACIAE.sh"             #
#                                                                              #
#                                               By An-Ke at CCNU on 10/17/2022 #
#                                                                              #
#                                                   Last updated on 02/14/2023 #
################################################################################
################################################################################


# Gets current date and time.
echo
current_date_and_time=$(date "+%Y%m%d%H%M%S")
echo "Now is $(date +"%Y-%m-%d %T")"
# mv ./usu.dat ./default_usu.dat
echo


################################################################################
################################################################################

# Modify the following statements as needed.

################################################################################
# For Makefile.
name_x="paciae.x"   # Name of the executable. User-defined.
name_f1="main_23"    # Names of the source files.
name_f2="parini_23"
name_f3="parcas_23"
name_f4="sfm_23"
name_f5="coales_23"
name_f6="hadcas_23"
name_f7="analy"
name_f8="p_23"
name_f9="stahad_23"
# name_f10="eps09"  # Not used now.
################################################################################


################################################################################
# Variable's alias for usu.dat

# Setups of the program-running.
# Total number of events in all program-runnings = n_run * n_eve.
n_run=1     # Number of program-running, i.e. total number of CPU cores used.
# n_run=${NP}   # Uncommented this line if submit jods on LSF system.
# ((n_run=SLURM_JOB_NUM_NODES*SLURM_NTASKS_PER_NODE))   # Uncommented for SLURM.
n_eve=10    # neve, number of events per program-running (per core).
n_out=1     # (D=neve/10) nout, outputs per n_out events.
n_osc=0     # (D=0) nosc, OSCAR-1997A/1999A/2013A output or not (1/2/3/0)
i_random=1  # (D=1) adj(26) selects random generator seed in PYTHIA.
            #         =0, default PYTHIA seed (19780503), can be used for debug
            #         =1, seed from the real-time system clock (h-min-s-ms)

# Setups of the incident particles.
na_proj=197 # nap, nucleon number of projectile
nz_proj=79  # nzp, proton number of projectile
na_targ=197 # nat, nucleon number of target
nz_targ=79  # nzt, proton number of target
            # For NN, NA (AN), and AB collisions
            #   p+p   : 1,1,1,1;
            #   p+pbar: 1,1,1,-1; pbar+p: 1,-1,1,1
            #   p+n   : 1,1,1, 0; n+p   : 1, 0,1,1
            #   n+n   : 1,0,1, 0;
            #   p+Pb  : 1,1,208,82; Pb+p: 208,82,1,1;
            #   Au+Au : 197,79,197,79; Pb+Pb: 208,82,208,82;
            #   Xe+Xe : 129,54,129,54; U + U: 238,92,238,92;
            #   Ag+Ag : 108,47,108,47; Cu+Cu: 63,29,63,29;
            #   Ru+Ru : 96,44,96,44; Zr+Zr: 96,40,96,40.

# Setups of the collision.
b_min=8.68  # bmin, min b parameter
b_max=9.87  # bmax, max b parameter
i_frame=1   # (D=1) collision frame, 1=collider, 0=fixed target.
energy=200  # (GeV) win, colliding CMS energy SQRT(s_NN) / incident momentum
i_channel=1 # (D=1) nchan, =0: inelastic (INEL)
            #               1: Non Single Difractive (NSD)
            #               2: qqb --> gamma^*/Z^0, used to generate Drell-Yan
            #               3: J/psi production
            #               4: heavy-flavor production
            #               5: direct photon
            #               6: soft only
            #               7: W+/- production
            #               8: pythia default (msel=1)
            #               9: Z0 production
            #               settings of 0,1,3,4,5,7,8 and 9 are ready
i_stage=4   # (D=4) adj(40), =1, stops event after parini
            #                 2, ... parcas
            #                 3, ... 
            #                 4, ... hadcas

# Setups of the parton initialization (parini).
k_pythia=1.5    # (D=1.5) adj(10), K factor multiplying the differential 
                #                    cross sections for hard parton-parton 
                #                    process in PYTHIA.

# Setups of the parton cascade/rescattering (parcas).
k_parcas=1      # (D=1) adj(1), K factor multiplying on the differential 
                #                 cross-section in parton cascade.
                #               =0, without parton cascade.

# Setups of the hadronization.
i_hadronization=0   # (D=0) adj(12), model for hadronization
                    #                =0, string fragmentation (sfm)
                    #                =1, Monte-Carlo coalescence (coal)
## String fragmentation parameter.
i_fragment_function=4   # (D=4 or 1) adj(29), selects fragmentation function in 
                        #                       sfm and coal.
                        #            For sfm, = MSTJ(11). Choice of longitudinal
                        #              fragmentation function, i.e. how large a 
                        #              fraction of the energy available a 
                        #              newly-created hadron takes.
                        #              =1: Lund symmetric fragmentation function
                        #              =2: Field–Feynman + Peterson/SLAC
                        #              =3: Lund + Peterson/SLAC
                        #              =4: Lund + Bowler, default in PYTHIA
                        #              =5: as = 4, but interpolate for c and b
                        #            For coal, sampling daughter parton energy 
                        #              fraction z taking from mother in 'funcz'.
                        #              =1: Lund string fragmentation function
                        #              =2: Field-Feynmman
                        #              =3: null
                        #              =4: null
                        #              =5: null
a_lund=0.3  # (D=0.3)  adj(6), alpha in the LUND string fragmentation function.
b_lund=0.58 # (D=0.58) adj(7), beta.
## Coalescence parameter.
# i_fragment_function=1 # (D=1)
i_phase_constraint=0    # (D=0) adj(21), phase-space constraint in coalescence.
                        #                =0, without
                        #                =1, with complete phase-space
                        #                =2, with position-space only
                        #                =3, with momentum-space only

# Setups of the hadron cascade/rescattering (hadcas).
i_hadcas=1  # (D=1) kjp21, = 0, without hadron cascade; 1, with
################################################################################


# The above setups are enogh to run program normall.


################################################################################
# The following statements are optional. Usually, they are good enough.

#  The initial geometry and effects.
b_samp=2    # (D=2) psno, selects b parameter sampling method for [b_min,b_max].
            #             =0, fixed b
            #              1, systematic sampling
            #              2, random sampling
i_overlap=0 # (D=0) adj(30), initial nuclon distribution.
            #                =0, without more requirements
            #                =1, distributes the participant nucleons in 
            #                   the overlapping region forcely
i_shadow=1  # (D=1) adj(5), choice of nuclear shadowing for partons.
            #               =0, without nuclear shadowing,
            #               =1, Wang's nuclear shadowing (PLB 527(2002)85).
i_cme=1     # (D=1) adj(23), choice of chiral magnetic effect(CME) for partons.
            #                =0, without CME-induced charge separation mechanism
            #                =1, with
parecc=0.       # (D=0) parecc, parameter converting initial spatial space 
                #               eccentricity to final momentum space ellipticity
pT_perturb=0.1  # (D=0.1) smadel, small perturbation of pT ellipse from circle 
                #                   when generates transverse momentum 
                #                   according to a Gaussian distribution.

# Processes in parcas.
i_inelastic=0   # (D=0) iparres, consider inelastic processes or not
                #                =0, elastic processes 1, 2, 3, 5, 8, and 9 only
                #                =1, elastic + inelastic
i_inel_proc=7   # (D=7) mstj1_1, when i_inelastic=1 (iparres=1),
                #                =6, with all the inelastic processes 4, 6, 7
                #                =7, only the inelastic process 7
i_t_shower=0    # (D=0) mstj1_2, when i_inelastic=1 (iparres=1),
                #                =0, without final state time-like parton shower
                #                =1, with

# Options in string fragmentation.
i_string_tension=4  # (D=4) kjp22, selects the form of string rension in sfm.
                    #              =1, variable single string tension
                    #               2, variable multiple string tension
                    #               3, variable (single+multiple) string tension
                    #               4, default constant string tension
cp0=1.              # (D=1.) cp0, cr0, parameters in parameterization 
cr0=0.2             # (D=0.2)            of multiple string effect
seco=0.05           # (D-0.05) seco, correction of PARJ(1) in popcorn mechanism
parj4=0.05          # (D=0.05) dparj4, PARJ(4) in PYTHIA

# Time accuracy.
dt_parini=0.00001   # (D=0.00001) ddt, min distinguishble collision time used in 
                    #                    the partonic initiation
dt_parcas=0.03      # (D=0.03) adj(19), time accuracy used in the parton cascade
dt_hadcas=0.1       # (D=0.1) adj(11), time accuracy used in the hadron cascade

# Statistical histogram bins in analy. The "asd(i)"" in usu.dat.
bin_1=0.35  # (D=0.35) y bin
bin_2=0.5   # (D=0.5)  pT bin
bin_3=0.35  # (D=0.35) eta bin
bin_4=0.3   # (D=0.3)  eT bin
bin_5=25    # (D=25)   multiplicity bin
i_y_or_eta=1    # (D=1) parp22, selects y or eta in partial-space statistics.
                #       = 0 , y
                #       = 1 , eta
################################################################################

################################################################################
# The following lines do not need to change. Normally automatical calculation.

# Total number of events in all program-runnings.
((tot_eve=n_run*n_eve))

# The "ipden" and "itden" in PACIAE for different systems.
i_proj=1
i_targ=1
if [[ "na_proj" = "1" ]]; then      # For N.
    i_proj=0
elif [[ "na_proj" > "1" ]]; then    # For A.
    i_proj=1
fi
if [[ "na_targ" = "1" ]]; then      # For N.
    i_targ=0
elif [[ "na_targ" > "1" ]]; then    # For A.
    i_targ=1
fi
# If others, please select them manully.
# i_proj=1  # ipden =0,  if projectile is nucleon (anti-nucleon)
            #       =1,  if nucleus
            #       =2,  for e+e-
            #       =11, if projectile is e- (e+)
            #       =12, if projectile is nu_e (nu_ebar)
            #       =13, if projectile is mu- (mu+)
            #       =14, if projectile is nu_mu (nu_mubar)
            #       =15, if projectile is tau- (tau+)
            #       =16, if projectile is nu_tau (nu_taubar)
# i_targ=1  # itden =0,  if target is nucleon (anti-nucleon)
            #       =1,  if nucleus
            #       =2,  for e+e-
            #...
            # for eA, nu_eA, etc.:
            # e^-A:     nap=1, nzp=-1, ipden=11, itden=1,
            # e^+A:     nap=1, nzp=1,  ipden=11, itden=1,
            # nu_eA:    nap=1, nzp=-1, ipden=12, itden=1,
            # nu_ebarA: nap=1, nzp=1,  ipden=12, itden=1,
################################################################################

# More options could be modified on the following "usu.dat" directly.
################################################################################
################################################################################
####################               usu.dat              ########################
echo "${n_eve},${n_out},${n_osc}                  ! neve,nout,nosc"   > usu.dat
echo "${na_proj},${nz_proj},${na_targ},${nz_targ} ! nap,nzp,nat,nzt" >> usu.dat
echo "${dt_parini},0.01,${b_min},${b_max},10      ! ddt,dtt,bmin,bmax,nmax"   >> usu.dat
echo "${i_hadcas},${i_frame},1.2,200.,1     ! kjp21,ifram,para7,para10,kjp20" >> usu.dat
echo "3.1416,${i_proj},${i_targ}                  ! pio,ipden,itden" >> usu.dat
echo "20,5,2                                 ! ispmax,isdmax,iflmax" >> usu.dat
echo "211,-211,321,-321,3122,-3122,3312,-3312,3334,-3334   ! KF code of particles" >> usu.dat
echo "2212,-2212,2112,-2112,3212,-3212,3112,3222,310,333   ! KF code of particles" >> usu.dat
echo "${bin_1},${bin_2},${bin_3},${bin_4},${bin_5}  ! asd(i), i=1,5" >> usu.dat
echo "-1,1                        ! afl(j=1,i=1,1), afl(j=1,i=1,2) " >> usu.dat
echo "0.,50.                      ! afl(j=1,i=2,1), afl(j=1,i=2,2) " >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "-1,1"                                                          >> usu.dat
echo "0.,50."                                                        >> usu.dat
echo "4.5,${i_y_or_eta},${energy}               ! parp21,parp22,win" >> usu.dat
echo "0.,0.5,1,1,${i_channel}       ! ttaup,taujp,iabsb,iabsm,nchan" >> usu.dat
echo "7.2,4.,${b_samp},40.,20.,0.,0.1 ! para13,para14,psno,para15,para16,ajpsi,vneum" >> usu.dat
echo "40.,40.,25.,10.                 ! para1_1,para1_2,para2,para4" >> usu.dat
echo "0.5,0.6,2.,0.05,0.35,20     ! tdh,cptl,cptu,cptl2,cptu2,itnum" >> usu.dat
echo "1,${i_inel_proc},${i_t_shower},2,3,2                     ! mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw" >> usu.dat
echo "${k_parcas},0.47,0.4,1000,${i_shadow},${a_lund},${b_lund},4,1.4,${k_pythia}   ! adj1(1)- adj1(10)  " >> usu.dat
echo "${dt_hadcas},${i_hadronization},26,18,1.,1,2.,0,${dt_parcas},1                ! adj1(11)- adj1(20) " >> usu.dat
echo "${i_phase_constraint},4.,${i_cme},0.15,0.4,${i_random},800000.,1.,${i_fragment_function},${i_overlap} ! adj1(21)- adj1(30) " >> usu.dat
echo "0.1,0.3,0.4,0.36,1.,0.,100.,3.,2.,${i_stage}  ! adj1(31)- adj1(40) " >> usu.dat
echo "${i_string_tension},2,2,0.,0      ! kjp22,kjp23,kjp24,parp78,mstptj" >> usu.dat
echo "${parecc},${i_inelastic},${pT_perturb},${parj4},${cp0},${cr0},${seco}   ! parecc,iparres,smadel,dparj4,cp0,cr0,seco" >> usu.dat
echo "0.35776,0.71552                 ! csp_31,csp_32"               >> usu.dat
echo "0.29198,0.58396,0.85464         ! csp_41,csp_42,csp_43"        >> usu.dat
echo "0.23894,0.47788,0.71113,0.91086 ! csp_51,csp_52,csp_53,csp_54" >> usu.dat
echo "0.19964,0.39927,0.59872,0.79707,0.99180 ! csp_61,csp_62,csp_63,csp_64,csp_65"   >> usu.dat
####################               usu.dat              ########################
################################################################################
#######                                                              ###########
#######   Annotation and KF code list are at the end of this file.   ###########
#######                                                              ###########
################################################################################



# The following statements are not required to be modified.



################################################################################
################################################################################
####################              Makefile              ########################
echo "# This is a toy Makefile for PACIAE."       > Makefile
echo "#     By Anke at CCNU on 10/17/2022 "      >> Makefile
echo "# How to use:"                             >> Makefile
echo "#   1. Type \"make\" command to compile and build PACIAE running file (${name_x})."   >> Makefile
echo "#   2. Type \"make clean\" command to clean the *.o , *.mod and *.x files. "          >> Makefile
echo "#   Usually, second command is not required to be used."      >> Makefile
echo                                                                >> Makefile
echo                                                                >> Makefile
echo "# The name of the executeble file."                           >> Makefile
# echo "target := paciae.x"                                         >> Makefile
echo "target := ${name_x}"                                          >> Makefile
echo                                                                >> Makefile
echo "# The names of source files."                                 >> Makefile
echo "#   Modify the following name to the needed one."             >> Makefile
echo "src_1 := ${name_f1}"                                          >> Makefile
echo "src_2 := ${name_f2}"                                          >> Makefile
echo "src_3 := ${name_f3}"                                          >> Makefile
echo "src_4 := ${name_f4}"                                          >> Makefile
echo "src_5 := ${name_f5}"                                          >> Makefile
echo "src_6 := ${name_f6}"                                          >> Makefile
echo "src_7 := ${name_f7}"                                          >> Makefile
echo "src_8 := ${name_f8}"                                          >> Makefile
echo "src_9 := ${name_f9}"                                          >> Makefile
echo "# scr_10 := ${name_f10}"                                      >> Makefile
echo                                                                >> Makefile
echo "# The intermediate files."                                    >> Makefile
echo "objects := \$(src_1).o \$(src_2).o \$(src_3).o \\"            >> Makefile
echo "           \$(src_4).o \$(src_5).o \$(src_6).o \\"            >> Makefile
echo "           \$(src_7).o \$(src_8).o \$(src_9).o"               >> Makefile
echo "#          \$(src_10).o"                                      >> Makefile
echo                                                                >> Makefile
echo "# Compiler"                                                   >> Makefile
echo "compiler := gfortran"                                         >> Makefile
echo "# Compiling flags"                                            >> Makefile
echo "comp_flags := -g -fbounds-check"                              >> Makefile
echo "# comp_flags := -g -Wall -fbounds-check"                      >> Makefile
echo "# comp_flags = -g -Wall -fbounds-Check -mcmodel=large"        >> Makefile
echo                                                                >> Makefile
echo "# Generating"                                                 >> Makefile
echo "\$(target) : \$(objects)"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -o \$(target) -O \$(objects)"   >> Makefile
echo                                                                >> Makefile
echo "\$(src_1).o: \$(src_1).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_1).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_2).o: \$(src_2).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_2).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_3).o: \$(src_3).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_3).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_4).o: \$(src_4).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_4).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_5).o: \$(src_5).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_5).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_6).o: \$(src_6).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_6).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_7).o: \$(src_7).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_7).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_8).o: \$(src_8).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_8).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_9).o: \$(src_9).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_9).f"              >> Makefile
echo                                                                >> Makefile
echo "# \$(src_10).o: \$(src_10).f"                                 >> Makefile
echo "# 	\$(compiler) \$(comp_flags) -c -O \$(src_10).f"         >> Makefile
echo                                                                >> Makefile
echo "# Cleans .o, .mod and .x files."                              >> Makefile
echo ".PHONY : clean"                                               >> Makefile
echo "clean : "                                                     >> Makefile
echo "	rm -rf *.o *.mod *.x"                                       >> Makefile
####################              Makefile              ########################
################################################################################
################################################################################



################################################################################
################################################################################
####################           Command lines            ########################
# Creates src folder.
if [ -d "./src/" ]; then
    echo "Folder src already exists."
else
    mkdir ./src
fi
f_file=$(ls *.f 2> /dev/null | wc -l)
if [ "${f_file}" != "0" ] ; then
    mv -f *.f ./src
fi

mv -f Makefile ./src
wait
# Enters src folder then compiles and builds executable file.
cd ./src 
echo
make
cd ../

# Creates bin folder.
if [ -d "./bin/" ]; then
    echo "Folder bin already exists."
else
    mkdir ./bin
fi
cp -f ./src/${name_x} ./bin

# Creates etc folder.
if [ -d "./etc/" ]; then
    echo "Folder etc already exists."
else
    mkdir ./etc
fi
# if [ -f "./EPS09*" ]; then
#    mv -f EPS09* ./etc
# fi
if [ -f "./default_usu.dat*" ]; then
   mv -f ./default_usu.dat ./etc
fi
mv -f usu.dat ./etc

# Creates log folder.
if [ -d "./log/" ]; then
    echo "Folder log already exists."
else
    mkdir ./log
fi
# Good for normal LINUX and SLURM, but not for LSF.
# TODO(Lei20221016): need to improve for LSF.
out_file=$(ls *.out 2> /dev/null | wc -l)
if [ "${out_file}" != "0" ] ; then
    mv -f *.out ./log
fi
log_file=$(ls *.log 2> /dev/null | wc -l)
if [ "${log_file}" != "0" ] ; then
    mv -f *.log ./log
fi
# mv -f *.out *.log ./log

# File name
if [[ "${na_proj}" = "1" && "${nz_proj}" = "1" ]]; then
    proj_name="p"
elif [[ "${na_proj}" = "1" && "${nz_proj}" = "-1" ]]; then
    proj_name="pbar"
elif [[ "${na_proj}" = "1" && "${nz_proj}" = "0" ]]; then
    proj_name="n"
elif [ "${na_proj}" = "63" ]; then
    proj_name="Cu${na_proj}"
elif [[ "${na_proj}" = "96" && "${nz_proj}" = "44" ]]; then
    proj_name="Ru${na_proj}"
elif [[ "${na_proj}" = "96" && "${nz_proj}" = "40" ]]; then
    proj_name="Zr${na_proj}"
elif [ "${na_proj}" = "108" ]; then
    proj_name="Ag${na_proj}"
elif [ "${na_proj}" = "129" ]; then
    proj_name="Xe${na_proj}"
elif [ "${na_proj}" = "197" ]; then
    proj_name="Au${na_proj}"
elif [ "${na_proj}" = "208" ]; then
    proj_name="Pb${na_proj}"
elif [ "${na_proj}" = "238" ]; then
    proj_name="U${na_proj}"
else
    proj_name="${na_proj}"
fi
# File name
if [[ "${na_targ}" = "1" && "${nz_targ}" = "1" ]]; then
    targ_name="p"
elif [[ "${na_targ}" = "1" && "${nz_targ}" = "-1" ]]; then
    targ_name="pbar"
elif [[ "${na_targ}" = "1" && "${nz_targ}" = "0" ]]; then
    targ_name="n"
elif [ "${na_targ}" = "63" ]; then
    targ_name="Cu${na_targ}"
elif [[ "${na_targ}" = "96" && "${nz_targ}" = "44" ]]; then
    proj_name="Ru${na_targ}"
elif [[ "${na_targ}" = "96" && "${nz_targ}" = "40" ]]; then
    proj_name="Zr${na_targ}"
elif [ "${na_targ}" = "108" ]; then
    targ_name="Ag${na_targ}"
elif [ "${na_targ}" = "129" ]; then
    targ_name="Xe${na_targ}"
elif [ "${na_targ}" = "197" ]; then
    targ_name="Au${na_targ}"
elif [ "${na_targ}" = "208" ]; then
    targ_name="Pb${na_targ}"
elif [ "${na_targ}" = "238" ]; then
    targ_name="U${na_targ}"
else
    targ_name="${na_targ}"
fi
# File name
if [ "${i_frame}" = "1" ]; then
    coll_name="COLL"
elif [ "${i_frame}" = "0" ]; then
    coll_name="FXT"
else
    coll_name="${i_frame}"
fi
# Creates files.
if [ -d "./${proj_name}_${targ_name}_${coll_name}/" ]; then
    echo "Folder ${proj_name}_${targ_name}_${coll_name} already exists."
else
    mkdir ./${proj_name}_${targ_name}_${coll_name}
fi
cd ./${proj_name}_${targ_name}_${coll_name}

if [ -d "./${energy}GeV/" ]; then
    echo "Folder ${energy}GeV already exists."
else
    mkdir ./${energy}GeV
fi
cd ./${energy}GeV
#ls -a

# b: impact parameter; 
# aL: a_lund
# bL: b_lund
# K : k_pythia (in parini)
#     kI: k_pythia (in parini)
#     kC: k_parcas
# iO: i_overlap
# iH: i_hadronization
# iHc: i_hadcas
# iStg: i_stage
# iC : i_channel
# iStr: i_string_tension
# iP : i_phase_constraint
# iFF : i_fragment_function
# iM : i_cMe
# pTp : pT_perturb
# Default 
dir="./b${b_min}_${b_max}_aL${a_lund}_bL${b_lund}_K${k_pythia}_iStg${i_stage}"
# dir="./b${b_min}_${b_max}_aL${a_lund}_bL${b_lund}_kI${k_pythia}_kC${k_parcas}_pTp${pT_perturb}_iO${i_overlap}_iH${i_hadronization}_iHc${i_hadcas}_iStr${i_string_tension}_iP${i_phase_constraint}_iFF${i_fragment_function}_iM${i_cme}_iC${i_channel}_iStg${i_stage}"

if [[ "${na_proj}" = "1" && "${na_targ}" = "1" ]]; then
# Elementary NN collision
    if [[ "${i_hadronization}" = "0" ]]; then
        dir="./Sfm_aL${a_lund}_bL${b_lund}_K${k_pythia}_iStg${i_stage}"
    elif [[ "${i_hadronization}" = "1" ]]; then
        dir="./Coal_K${k_pythia}_iStg${i_stage}"
    # else
    fi
# elif [[  ]]; then
else
# pA, AA collisions
    if [[ "${i_hadronization}" = "0" ]]; then
        dir="./b${b_min}_${b_max}_sfm_aL${a_lund}_bL${b_lund}_K${k_pythia}_iStg${i_stage}"
    elif [[ "${i_hadronization}" = "1" ]]; then
        dir="./b${b_min}_${b_max}_coal_K${k_pythia}_iStg${i_stage}"
    # else
    fi
fi


if [ -d "${dir}/" ]; then
    echo "Folder ${dir} already exists."
else
    mkdir ${dir}
fi
#ls -a

if [ -d "${dir}/${tot_eve}_events/" ]; then
    echo "Folder ${tot_eve}_events already exists."
else
    mkdir ${dir}/${tot_eve}_events
fi
cd ${dir}/${tot_eve}_events
echo
pwd
echo

# Not required to be modified. Pre-setups for CPU/program-running.
i_cpu=0      # Do not change i_cpu=0. ID od CPU (from 0)
i_run=1      # DO not change i_run=1. Counter for program-running.

while [ ${i_run} -le ${n_run} ]
do
    if [ -d "./PACIAE_${i_run}/" ]; then
        echo
    else
        mkdir ./PACIAE_${i_run}
    fi
    # ls -a
    cp -f ../../../../bin/${name_x} ./PACIAE_${i_run}/${i_run}_${name_x}
    cp -f ../../../../etc/usu.dat ./PACIAE_${i_run}/
    # cp -f ../../../../etc/EPS09* ./PACIAE_${i_run}/
    cd ./PACIAE_${i_run}/
    # ls -a
    nohup time ./${i_run}_${name_x} > screen.log &
# Binds one running task to one CPU.
    # nohup time taskset -c ${i_cpu} ./${i_run}_${name_x} > screen.log &
    cd ..
    # ls -a

    echo "PACIAE-${i_run}"

    ((i_cpu=i_cpu+1))
    ((i_run=i_run+1))
done
((i_run=i_run-1))

wait

cd ../../../..
echo
ls -a
echo
echo "PACIAE-script running ended!"
echo
####################           Command lines            ########################
################################################################################
################################################################################



################################################################################
################################################################################
######################        Annotation of usu.dat         ####################
# neve,nout,nosc (D=10, 1, 0)
# neve: events number to be generated
# nout: output the event per nout events
# nosc: OSCAR standard output (oscar.out)
#       = 0 : no OSCAR output, see subroutine oscar
#       = 1 : OSCAR1997A (final_id_p_x, just PACIAE final output)
#       = 2 : OSCAR1999A (full_event_history)
#       = 3 : OSCAR2013A (full_event_history, dummy now)
#
# nap,nzp,nat,nzt
# for NN, NA(AN), AA, etc.
#  nap(nzp): nucleons (protons) number of projectile 
#  nat(nzt): nucleons (protons) number of target 
#            for NN along with ipden=itden=0
#             p+p   : 1, 1, 1, 1 ;
#             p+pbar: 1, 1, 1,-1 ;
#             pbar+p: 1,-1, 1, 1 ;
#             p+n   : 1, 1, 1, 0 ;
#             n+p   : 1, 0, 1, 1 ;
#             n+n   : 1, 0, 1, 0 ;
#            for pA(Ap) along with ipden=0, itden=1 (ipden=1, itden=0)
#             p+Pb  : 1, 1, 208, 82;
#             Pb+p  : 208, 82, 1, 1;
#            for A+B along with ipden=itden=1 
#             Au+Au: 197, 79, 197, 79; Pb+Pb: 208, 82, 208, 82;
#             Xe+Xe: 129, 54, 129, 54; U + U: 238, 92, 238, 92;
#             Ag+Ag: 108, 47, 108, 47; Cu+Cu:  63, 29,  63, 29;
#             Ru+Ru:  96, 44,  96, 44; Zr+Zr:  96, 40,  96, 40.
# for eA, nu_eA, etc.
#  e^-A:     nap=1, nzp=-1, ipden=11, itden=1, 
#  e^+A:     nap=1, nzp= 1, ipden=11, itden=1, 
#  nu_eA:    nap=1, nzp=-1, ipden=12, itden=1,
#  nu_ebarA: nap=1, nzp= 1, ipden=12, itden=1.
#
# ddt,dtt,bmin,bmax,nmax (D=0.00001, 0.01, xxx, xxx, 10)
#  ddt: minimum distinguishble collision time interval used in 
#       partonic initiation in parini.f
#  dtt: not used now
#  bmin: minimum impact parameters, 
#  bmax: maximum impact parameters,
#  nmax: the number of intervals segmented in [bmin,bmax] when psno=1
#
# kjp21,ifram,para7,para10,kjp20 (D=1, 1, 1.2, 200, 1)
#  kjp21: =0, without hadron rescattering
#         =1, with hadron rescattering
#  ifram: choice collision system type
#         =0, fixed target
#         =1, collider
#  para7: proper formation time in rest-frame of particle
#  para10: largest allowed size of partonic (hadronic) rescattering
#          region which is product of para10 and target radius
#  kjp20: choice the cross sections in hadron rescattering (hadcas.f)
#         =1, constant cross sections 
#         =0, energy dependent cross sections
#
# pio,ipden,itden (D=3.1416, 1, 1 for A+B)
#  pio: pi=3.1416
#  ipden: =0, if projectile is nucleon (anti-nucleon)
#         =1, if projectile is nucleus
#         =2, for e+e-
#         =11, if projectile is e- (e+)
#         =12, if projectile is nu_e (nu_ebar)
#         =13, if projectile is mu- (mu+)
#         =14, if projectile is nu_mu (nu_mubar)
#         =15, if projectile is tau- (tau+)
#         =16, if projectile is nu_tau (nu_taubar)
#  itden: =0, if target is nucleon (anti-nucleon)
#         =1, if projectile is a nucleus
#         =2, for e+e-
#         ...
#
# ispmax,isdmax,iflmax   (D=20, 5, 2)
#  ispmax: maximum # of different particle pieces to be considered
#  isdmax: maximum # of different distributions to be considered
#  iflmax: maximum # of windows to be set, =0 means no window at all
#
# ispkf(i,i=1,ispmax):
# KF code: particle code used in PYTHIA and PACIAE, 
#          (list at the end and see detail in reference: arXiv:hep-ph/0603175)
#
# asd(i=1,isdmax): interval of the i-th distribution
#  for pp, pbarp, pA(Ap), AB etc. 
#      (D=0.35, 0.5, 0.35, 0.3, 25)
#      i=1: rapidity distribution
#       =2: transverse monmentum distribution
#       =3: pesudorapidity distribution
#  for ep, nu_ep, etc. 
#      i=1: Q^2=-q^2 (fq2 in code) distribution
#       =2: W^2 (w21) distribution
#       =3: y (yyl) distribution
#       =4: p_h (pph) distribution
#       =5: z (zl) distribution
#
# afl(j,i,1): lower-boundary of i-th window for j-th particle
# afl(j,i,2): upper-boundary  of i-th window for j-th particle
#  for pp, pbarp, pA(Ap), AB etc. 
#      i=1, rapidity/pesudorapidity window (D= -1.,1. )
#       =2, transverse monmentum           (D= 0.,50. )
#  for ep, nu_ep, etc. 
#      i=1, Q^2=-q^2 window
#       =2, W^2
#       =3, y
#       =4, p_h (haron momentum)
#       =5: z
#
# parp21,parp22,win (D=4.5, 10., xxx)
#  parp21: lowest CM energy running 'pythia'
#  parp22: select y or eta in partial phase-space statistics (analy.f)
#          = 0 , y
#          = 1 , eta
#  win = cms energy if ifram=1 (collider)
#      = incident momentum if ifram=0 (fixed target)
#
# ttaup,taujp,iabsb,iabsm,nchan (D=0., 0.5, 1, 1, 1)
#  ttaup: proper formation time of particles generated in hadronic rescattering
#  taujp: formation time of J/Psi
#  iabsb: =0, without J/Psi (Psi') + baryon
#         =1, with J/Psi (Psi') + baryon
#  iabsm: =0, without J/Psi (Psi') + meson
#         =1, with J/Psi (Psi') + meson
#  nchan: to choose which subset of parton-parton subprocesses to include in
#         the generration
#         =0, inelastic (INEL)
#         =1, Non Single Difractive (NSD) 
#         =2, qqb --> gamma^*/Z^0, to generate Drell-Yan
#         =3, J/psi production
#         =4, heavy-flavor production
#         =5, direct photon
#         =6, soft only
#         =7, default PYTHIA 
#
# para13,para14,psno,para15,para16,ajpsi,vneum (D=7.2, 4., 2, 40., 20., 0, 0.1)
#  para13: totle cross-section of J/Psi + n
#  para14: totle cross-section of J/Psi + meson 
#  psno: =0 fixed impact parameter 
#        =1 impact parameter is sampled by systematic sampling method
#        =2 randomly sampled impact parameter 
#  para15: totle cross-section of Psi' + n
#  para16: totle cross-section of Psi' + meson
#  ajpsi: not used now
#  vneum: relevant to average binary collision number, now it is recalculated
#         in program
#
# para1_1,para1_2,para2,para4 (D= 40, 40, 25, 10 for Au+Au; 75, 75... for Pb+Pb)
#  para1_1: nn total cross section used in parton initiation
#  para1_2: nn total cross section used in hadron cascade
#  para2: totle cross-section of pi-nucleon
#  para4: totle cross-section of pi-pi
#
# tdh,cptl,cptu,cptl2,cptu2,itnum (D=0.5, 0.6, 2., 0.05, 0.35, 20)
#  tdh: not used now
#  cptl: not used now
#  cptu: not used now
#  cptl2: not used now
#  cptu2: not used now
#  itnum: not used now
#
# mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw (D=1, 7, 0, 2, 3, 2)
#  mstu21: parameter mstu(21) in PYTHIA
#  mstj1_1: = 6, with inelastic processes 4, 6, and 7 if iparres=1 (in parcas.f)
#           = 7, with inelastic process 7 only if iparres=1 (in parcas.f)
#  mstj1_2: = 0, w/o final state time-like parton shower if iparres=1
#           = 1, w/ final state time-like parton shower if iparres=1
#  mstj2: parameter mstj(2) in PYTHIA
#  mstj3: parameter mstj(3) in PYTHIA
#  itorw: =1 executing pyevnt, =2 executing pyevnw
#
# adj1(i), i=1,40: switches and/or parameters
# ------------------------------------------------------------------------------
#     D= 1-10 :  1., 0.47, 0.4, 1000,   1, 0.3,    0.58,  4,  1.4, 1.5
#        11-20: 0.1,    0,  26,   18,  1.,   1,     2 .,  0, 0.03, 1
#        21-30:   0,   4.,   1, 0.15, 0.4,   1, 800000., 1.,    4, 0
#        31-40: 0.1,  0.3, 0.4, 0.36,  1.,  0.,    100., 3.,   2., 4.
# ------------------------------------------------------------------------------
#     i=1: K factor in parton cascade
#       2: $\alpha_s$, effective coupling constant in parton rescattering
#       3: parameter \mu^2 (tcut in program), the regulation factor introduced 
#          in the parton-parton differential cross section   
#       4: parameter `idw', the number of intervals in the numerical integration
#       5: =1: with Wang's nuclear shadowing (PLB 527(2002)85),
#          =0: without nuclear shadowing
#       6: \alpha in the LUND string fragmentation function (parj(41) in PYTHIA)
#       7: $\beta$ in the LUND string fragmentation function ( parj(42) )
#       8: mstp(82) in PYTHIA
#       9: parp(81) in PYTHIA
#       10: K factor (parp(31) in PYTHIA).
#       11: time accuracy used in the hadron cascade
#       12: model for hadronization 
#           =0, string fragmentation,
#           =1, Monte Carlo coalescence model,
#           =2, statistical hadronization. 
#       13: dimension of meson table considered in coalescence model.
#       14: dimension of baryon table considered coalescence model.
#       15: string tension of qqbar simple string
#       16: number of loops in the deexcitation of energetic quark in the 
#           Monte Carlo coalescence model.
#       17: the threshold energy in the deexcitation of energetic quark in 
#           the Monte Carlo coalescence model.
#       18: =0, rest partons hadronize by string fragmentation
#           =1, rest partons hadronize by coalescence
#       19: time accuracy used in the parton cascade 
#       20: the optional parton-parton cross section in the parton rescattering
#           =0, LO pQCD parton-parton cross section 
#           =1, keeping only leading divergent terms in the LO pQCD 
#               parton-parton cross section (B. Zhang)
#           =2, the same as 0 but flat scattering angle distribution is assumed
#           =3, the same as 1 but flat scattering angle distribution is assumed.
#       21: with or without phase space constraint in the MC coalescence model
#           =0, without phase space constraint
#           =1, with complete phase space constraint
#           =2, with spatial phase space constraint
#           =3, with momentum phase space constraint
#       22: critical value (D=4) of the product of radii in position and 
#           momentum phase spaces 
#       23: switch for chiral magnetic effect (CME)
#           =0: CME off
#           =1: CME on
#       24: the virtuality cut ('tl0' in program) in the time-like radiation in  
#           parton rescattering 
#       25: $\Lambda_{QCD}$ in PYTHIA 
#       26: selection of random number seed,
#           =0, default PYTHIA seed (19780503), can be used for debug
#           =other, seed from the real-time colock
#       27: largest momentum allowed for produced particle
#       28: concerned to the largest position allowed for produced particle
#       29: For sfm in PYTHIA, it is MSTJ(11). Choice of longitudinal 
#            fragmentation function, i.e. how large a fraction of the energy 
#            available a newly-created hadron takes.
#           =1: Lund symmetric fragmentation function, see PARJ(41) - PARJ(45)
#           =2: Field–Feynman + Peterson/SLAC, see PARJ(51) PARJ(59)
#           =3: Lund + Peterson/SLAC (light flavor + heavier)
#           =4: default PYTHIA. Lund + Bowler.
#           =5: as = 4, but interpolate for c and b; see PARJ(46) and PARJ(47).
#           For coal, sampling daughter parton energy fraction z taking from 
#            mother in 'funcz' for coal.
#           =1: by Lund string fragmentation function
#           =2: by Field-Feynmman fragmentation function
#           =3: null
#           =4: null
#           =5: null
#       30: =1, distribute the participant nucleons in overlapping region
#           =0, without more requirements
#       31: parj(1) in PYTHIA 
#       32: parj(2) in PYTHIA 
#       33: parj(3) in PYTHIA 
#       34: parj(21) in PYTHIA 
#       35: mstp(91) in PYTHIA, parton transverse momentum ($k_{\perp}$)  
#           distribution inside hadron
#           =1: Gaussian
#           =2: exponential
#       36: w/ or w/o phenomenological parton energy loss in parton rescattering
#           =0, without 
#           =1, with 
#       37: the coefficient in phenomenological parton energy loss
#       38: $p_T$ cut in phenomenological parton energy loss
#       39: width of Gaussian parton k_perp distribution in hadron if mstp(91)=1
#           width of exponential k_perp distribution in hadron if mstp(91)=2
#       40: optional event stopping point 
#          =1, after parton initiation,
#          =2, after parton rescattering,
#          =4, after hadron rescattering 
#
# kjp22,kjp23,kjp24,parp78,mstptj (D=4, 2, 2, 0., 0)
#  kjp22: =1, variable single string tension and parj(1) etc.
#         =2, variable multiple string tension and parj(1) etc.
#         =3, variable (single+multiple) string tension and parj(1) etc.
#         =4, default string tension and parj(1) etc.
#  kjp23: optional model for the calculation of participant nucleon # (npart)
#         =1, geometric model
#         =2, Glauber model
#  kjp24: optional distribution in Glauber model
#         =1, sharp sphere
#         =2, Woods-Saxon
#  parp78: parameter controling amount of colour reconnection in final state
#  mstptj: =0, input mstp(111) (mstj(1)) for pp, pA (AP), and AA (for e+e-) in
#              PACIAE simulation devaloped from partonic initial stage, 
#              to partonic rescattering, hadronization, and to hadronic
#              rescttering stage  
#          =1, PYTHIA like simulation without partonic & hadronic 
#              rescatterings but with setting of kjp21=0
#
# parecc,iparres,smadel,dparj4,cp0,cr0,seco (D=0., 0, 0.1, 0.05, 1., 0.2, 0.05)
#  parecc: proportional factor between initial spatial space eccentricity and 
#          final momentum space ellipticity 
#  iparres: =0 consider elastic parton-parton collisions only in parton cascade
#           =1 otherwise
#  smadel: small perpurbation of ellipse from circle
#  dparj4: default parj(4)
#  cp0,cr0: parameters in parameterization of multiple string effect
#  seco: parameter in popcorn mechanism for correction of parj(1)
#
#  cumulent sum of flavor generation probability, for an example: 
#  p_31=(1-amd/sms)/2., p_32=(1-amu/sms)/2., p_33=(1-ams/sms)/2.; 
#  amd (amu, ams): d (u,s) quark constituent mass, sms=amd+amu+ams; 
#  p_31+p_32+p33=1;
#  cumulent sum: csp_31=p_31, csp_32=p_31+p_32
#  u and d have same constituent mass and probability
######################        Annotation of usu.dat         ####################
################################################################################



################################################################################
################################################################################
#
#                     List of KF codes in program
#
#             1    d                            -1    dbar            
#             2    u                            -2    ubar            
#             3    s                            -3    sbar            
#             4    c                            -4    cbar            
#             5    b                            -5    bbar            
#             6    t                            -6    tbar            
#             7    b'                           -7    b'bar           
#             8    t'                           -8    t'bar           
#            11    e-                          -11    e+              
#            12    nu_e                        -12    nu_ebar         
#            13    mu-                         -13    mu+             
#            14    nu_mu                       -14    nu_mubar        
#            15    tau-                        -15    tau+            
#            16    nu_tau                      -16    nu_taubar       
#            17    tau'-                       -17    tau'+           
#            18    nu'_tau                     -18    nu'_taubar      
#            21    g               
#            22    gamma           
#            23    Z0              
#            24    W+                          -24    W-              
#            25    h0              
#            32    Z'0             
#            33    Z"0             
#            34    W'+                         -34    W'-             
#            35    H0              
#            36    A0              
#            37    H+                          -37    H-              
#            39    Graviton        
#            41    R0                          -41    Rbar0           
#            42    LQ_ue                       -42    LQ_uebar        
#          2101    ud_0                      -2101    ud_0bar         
#          3101    sd_0                      -3101    sd_0bar         
#          3201    su_0                      -3201    su_0bar         
#          4101    cd_0                      -4101    cd_0bar         
#          4201    cu_0                      -4201    cu_0bar         
#          4301    cs_0                      -4301    cs_0bar         
#          5101    bd_0                      -5101    bd_0bar         
#          5201    bu_0                      -5201    bu_0bar         
#          5301    bs_0                      -5301    bs_0bar         
#          5401    bc_0                      -5401    bc_0bar         
#          1103    dd_1                      -1103    dd_1bar         
#          2103    ud_1                      -2103    ud_1bar         
#          2203    uu_1                      -2203    uu_1bar         
#          3103    sd_1                      -3103    sd_1bar         
#          3203    su_1                      -3203    su_1bar         
#          3303    ss_1                      -3303    ss_1bar         
#          4103    cd_1                      -4103    cd_1bar         
#          4203    cu_1                      -4203    cu_1bar         
#          4303    cs_1                      -4303    cs_1bar         
#          4403    cc_1                      -4403    cc_1bar         
#          5103    bd_1                      -5103    bd_1bar         
#          5203    bu_1                      -5203    bu_1bar         
#          5303    bs_1                      -5303    bs_1bar         
#          5403    bc_1                      -5403    bc_1bar         
#          5503    bb_1                      -5503    bb_1bar         
#           111    pi0             
#           211    pi+                        -211    pi-             
#           221    eta             
#           311    K0                         -311    Kbar0           
#           130    K_L0            
#           310    K_S0            
#           321    K+                         -321    K-              
#           331    eta'            
#           411    D+                         -411    D-              
#           421    D0                         -421    Dbar0           
#           431    D_s+                       -431    D_s-            
#           441    eta_c           
#           511    B0                         -511    Bbar0           
#           521    B+                         -521    B-              
#           531    B_s0                       -531    B_sbar0         
#           541    B_c+                       -541    B_c-            
#           551    eta_b           
#           113    rho0            
#           213    rho+                       -213    rho-            
#           223    omega           
#           313    K*0                        -313    K*bar0          
#           323    K*+                        -323    K*-             
#           333    phi             
#           413    D*+                        -413    D*-             
#           423    D*0                        -423    D*bar0          
#           433    D*_s+                      -433    D*_s-           
#           443    J/psi           
#           513    B*0                        -513    B*bar0          
#           523    B*+                        -523    B*-             
#           533    B*_s0                      -533    B*_sbar0        
#           543    B*_c+                      -543    B*_c-           
#           553    Upsilon         
#         10113    b_10            
#         10213    b_1+                     -10213    b_1-            
#         10223    h_1             
#         10313    K_10                     -10313    K_1bar0         
#         10323    K_1+                     -10323    K_1-            
#         10333    h'_1            
#         10413    D_1+                     -10413    D_1-            
#         10423    D_10                     -10423    D_1bar0         
#         10433    D_1s+                    -10433    D_1s-           
#         10443    h_1c            
#         10513    B_10                     -10513    B_1bar0         
#         10523    B_1+                     -10523    B_1-            
#         10533    B_1s0                    -10533    B_1sbar0        
#         10543    B_1c+                    -10543    B_1c-           
#         10553    h_1b            
#         10111    a_00            
#         10211    a_0+                     -10211    a_0-            
#         10221    f_0             
#         10311    K*_00                    -10311    K*_0bar0        
#         10321    K*_0+                    -10321    K*_0-           
#         10331    f'_0            
#         10411    D*_0+                    -10411    D*_0-           
#         10421    D*_00                    -10421    D*_0bar0        
#         10431    D*_0s+                   -10431    D*_0s-          
#         10441    chi_0c          
#         10511    B*_00                    -10511    B*_0bar0        
#         10521    B*_0+                    -10521    B*_0-           
#         10531    B*_0s0                   -10531    B*_0sbar0       
#         10541    B*_0c+                   -10541    B*_0c-          
#         10551    chi_0b          
#         20113    a_10            
#         20213    a_1+                     -20213    a_1-            
#         20223    f_1             
#         20313    K*_10                    -20313    K*_1bar0        
#         20323    K*_1+                    -20323    K*_1-           
#         20333    f'_1            
#         20413    D*_1+                    -20413    D*_1-           
#         20423    D*_10                    -20423    D*_1bar0        
#         20433    D*_1s+                   -20433    D*_1s-          
#         20443    chi_1c          
#         20513    B*_10                    -20513    B*_1bar0        
#         20523    B*_1+                    -20523    B*_1-           
#         20533    B*_1s0                   -20533    B*_1sbar0       
#         20543    B*_1c+                   -20543    B*_1c-          
#         20553    chi_1b          
#           115    a_20            
#           215    a_2+                       -215    a_2-            
#           225    f_2             
#           315    K*_20                      -315    K*_2bar0        
#           325    K*_2+                      -325    K*_2-           
#           335    f'_2            
#           415    D*_2+                      -415    D*_2-           
#           425    D*_20                      -425    D*_2bar0        
#           435    D*_2s+                     -435    D*_2s-          
#           445    chi_2c          
#           515    B*_20                      -515    B*_2bar0        
#           525    B*_2+                      -525    B*_2-           
#           535    B*_2s0                     -535    B*_2sbar0       
#           545    B*_2c+                     -545    B*_2c-          
#           555    chi_2b          
#        100443    psi'            
#        100553    Upsilon'        
#          3122    Lambda0                   -3122    Lambdabar0      
#          4122    Lambda_c+                 -4122    Lambda_cbar-    
#          4132    Xi_c0                     -4132    Xi_cbar0        
#          4232    Xi_c+                     -4232    Xi_cbar-        
#          5122    Lambda_b0                 -5122    Lambda_bbar0    
#          5132    Xi_b-                     -5132    Xi_bbar+        
#          5232    Xi_b0                     -5232    Xi_bbar0        
#          5142    Xi_bc0                    -5142    Xi_bcbar0       
#          5242    Xi_bc+                    -5242    Xi_bcbar-       
#          5342    Omega_bc0                 -5342    Omega_bcbar0    
#          2112    n0                        -2112    nbar0           
#          2212    p+                        -2212    pbar-           
#          3112    Sigma-                    -3112    Sigmabar+       
#          3212    Sigma0                    -3212    Sigmabar0       
#          3222    Sigma+                    -3222    Sigmabar-       
#          3312    Xi-                       -3312    Xibar+          
#          3322    Xi0                       -3322    Xibar0          
#          4112    Sigma_c0                  -4112    Sigma_cbar0     
#          4212    Sigma_c+                  -4212    Sigma_cbar-     
#          4222    Sigma_c++                 -4222    Sigma_cbar--    
#          4312    Xi'_c0                    -4312    Xi'_cbar0       
#          4322    Xi'_c+                    -4322    Xi'_cbar-       
#          4332    Omega_c0                  -4332    Omega_cbar0     
#          4412    Xi_cc+                    -4412    Xi_ccbar-       
#          4422    Xi_cc++                   -4422    Xi_ccbar--      
#          4432    Omega_cc+                 -4432    Omega_ccbar-    
#          5112    Sigma_b-                  -5112    Sigma_bbar+     
#          5212    Sigma_b0                  -5212    Sigma_bbar0     
#          5222    Sigma_b+                  -5222    Sigma_bbar-     
#          5312    Xi'_b-                    -5312    Xi'_bbar+       
#          5322    Xi'_b0                    -5322    Xi'_bbar0       
#          5332    Omega_b-                  -5332    Omega_bbar+     
#          5412    Xi'_bc0                   -5412    Xi'_bcbar0      
#          5422    Xi'_bc+                   -5422    Xi'_bcbar-      
#          5432    Omega'_bc0                -5432    Omega'_bcba     
#          5442    Omega_bcc+                -5442    Omega_bccbar-   
#          5512    Xi_bb-                    -5512    Xi_bbbar+       
#          5522    Xi_bb0                    -5522    Xi_bbbar0       
#          5532    Omega_bb-                 -5532    Omega_bbbar+    
#          5542    Omega_bbc0                -5542    Omega_bbcbar0   
#          1114    Delta-                    -1114    Deltabar+       
#          2114    Delta0                    -2114    Deltabar0       
#          2214    Delta+                    -2214    Deltabar-       
#          2224    Delta++                   -2224    Deltabar--      
#          3114    Sigma*-                   -3114    Sigma*bar+      
#          3214    Sigma*0                   -3214    Sigma*bar0      
#          3224    Sigma*+                   -3224    Sigma*bar-      
#          3314    Xi*-                      -3314    Xi*bar+         
#          3324    Xi*0                      -3324    Xi*bar0         
#          3334    Omega-                    -3334    Omegabar+       
#          4114    Sigma*_c0                 -4114    Sigma*_cbar0    
#          4214    Sigma*_c+                 -4214    Sigma*_cbar-    
#          4224    Sigma*_c++                -4224    Sigma*_cbar--   
#          4314    Xi*_c0                    -4314    Xi*_cbar0       
#          4324    Xi*_c+                    -4324    Xi*_cbar-       
#          4334    Omega*_c0                 -4334    Omega*_cbar0    
#          4414    Xi*_cc+                   -4414    Xi*_ccbar-      
#          4424    Xi*_cc++                  -4424    Xi*_ccbar--     
#          4434    Omega*_cc+                -4434    Omega*_ccbar-   
#          4444    Omega*_ccc++              -4444    Omega*_cccbar-  
#          5114    Sigma*_b-                 -5114    Sigma*_bbar+    
#          5214    Sigma*_b0                 -5214    Sigma*_bbar0    
#          5224    Sigma*_b+                 -5224    Sigma*_bbar-    
#          5314    Xi*_b-                    -5314    Xi*_bbar+       
#          5324    Xi*_b0                    -5324    Xi*_bbar0       
#          5334    Omega*_b-                 -5334    Omega*_bbar+    
#          5414    Xi*_bc0                   -5414    Xi*_bcbar0      
#          5424    Xi*_bc+                   -5424    Xi*_bcbar-      
#          5434    Omega*_bc0                -5434    Omega*_bcbar0   
#          5444    Omega*_bcc+               -5444    Omega*_bccbar-  
#          5514    Xi*_bb-                   -5514    Xi*_bbbar+      
#          5524    Xi*_bb0                   -5524    Xi*_bbbar0      
#          5534    Omega*_bb-                -5534    Omega*_bbbar+   
#          5544    Omega*_bbc0               -5544    Omega*_bbcbar0  
#          5554    Omega*_bbb-               -5554    Omega*_bbbbar+  
#       1000001    ~d_L                   -1000001    ~d_Lbar         
#       1000002    ~u_L                   -1000002    ~u_Lbar         
#       1000003    ~s_L                   -1000003    ~s_Lbar         
#       1000004    ~c_L                   -1000004    ~c_Lbar         
#       1000005    ~b_1                   -1000005    ~b_1bar         
#       1000006    ~t_1                   -1000006    ~t_1bar         
#       1000011    ~e_L-                  -1000011    ~e_L+           
#       1000012    ~nu_eL                 -1000012    ~nu_eLbar       
#       1000013    ~mu_L-                 -1000013    ~mu_L+          
#       1000014    ~nu_muL                -1000014    ~nu_muLbar      
#       1000015    ~tau_1-                -1000015    ~tau_1+         
#       1000016    ~nu_tauL               -1000016    ~nu_tauLbar     
#       1000021    ~g              
#       1000022    ~chi_10         
#       1000023    ~chi_20         
#       1000024    ~chi_1+                -1000024    ~chi_1-         
#       1000025    ~chi_30         
#       1000035    ~chi_40         
#       1000037    ~chi_2+                -1000037    ~chi_2-         
#       1000039    ~Gravitino      
#       2000001    ~d_R                   -2000001    ~d_Rbar         
#       2000002    ~u_R                   -2000002    ~u_Rbar         
#       2000003    ~s_R                   -2000003    ~s_Rbar         
#       2000004    ~c_R                   -2000004    ~c_Rbar         
#       2000005    ~b_2                   -2000005    ~b_2bar         
#       2000006    ~t_2                   -2000006    ~t_2bar         
#       2000011    ~e_R-                  -2000011    ~e_R+           
#       2000012    ~nu_eR                 -2000012    ~nu_eRbar       
#       2000013    ~mu_R-                 -2000013    ~mu_R+          
#       2000014    ~nu_muR                -2000014    ~nu_muRbar      
#       2000015    ~tau_2-                -2000015    ~tau_2+         
#       2000016    ~nu_tauR               -2000016    ~nu_tauRbar     
#       3000111    pi_tc0          
#       3000211    pi_tc+                 -3000211    pi_tc-          
#       3000221    pi'_tc0         
#       3000331    eta_tc0         
#       3000113    rho_tc0         
#       3000213    rho_tc+                -3000213    rho_tc-         
#       3000223    omega_tc        
#       3100021    V8_tc           
#       3100111    pi_22_1_tc      
#       3200111    pi_22_8_tc      
#       3100113    rho_11_tc       
#       3200113    rho_12_tc       
#       3300113    rho_21_tc       
#       3400113    rho_22_tc       
#       4000001    d*                     -4000001    d*bar           
#       4000002    u*                     -4000002    u*bar           
#       4000011    e*-                    -4000011    e*bar+          
#       4000012    nu*_e0                 -4000012    nu*_ebar0       
#       5000039    Graviton*       
#       9900012    nu_Re           
#       9900014    nu_Rmu          
#       9900016    nu_Rtau         
#       9900023    Z_R0            
#       9900024    W_R+                   -9900024    W_R-            
#       9900041    H_L++                  -9900041    H_L--           
#       9900042    H_R++                  -9900042    H_R--           
#       9900110    rho_diff0       
#       9900210    pi_diffr+              -9900210    pi_diffr-       
#       9900220    omega_di        
#       9900330    phi_diff        
#       9900440    J/psi_di        
#       9902110    n_diffr0               -9902110    n_diffrbar0     
#       9902210    p_diffr+               -9902210    p_diffrbar-     
#       9900443    cc~[3S18]       
#       9900441    cc~[1S08]       
#       9910441    cc~[3P08]       
#       9900553    bb~[3S18]       
#       9900551    bb~[1S08]       
#       9910551    bb~[3P08]       
################################################################################
################################################################################
################################################################################
#                                                                              #
#                                                                              #
################################################################################
################################################################################
#                                                                              #
#         PPPPPPP       AAAAA      CCCCCCC     IIIIIII     AAAAA     EEEEEEEEE #
#        P      p     A     A    C       C       I       A     A    E          #
#       P       p   A       A   C               I      A       A   E           #
#      P      p    A       A   C               I      A       A   E            #
#     PPPPPPP     AAAAAAAAA   C               I      AAAAAAAAA   EEEEEEEEE     #
#    P           A       A   C               I      A       A   E              #
#   P           A       A   C               I      A       A   E               #
#  P           A       A    C       C      I      A       A   E                #
# P           A       A     CCCCCCC    IIIIIII   A       A   EEEEEEEEE         #
#                                                                              #
################################################################################
################################################################################
