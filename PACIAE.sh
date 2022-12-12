#!/bin/bash
#SBATCH -J PACIAE_SLURM       # Task running name. Arbitrary.
#SBATCH -N 1                  # Number of nodes. Modify as needed.
#SBATCH --ntasks-per-node=1   # Number of cores used per node. Modify as needed.
#SBATCH -p middle             # Name of sequence. Modify as needed.

# Above statements are for normal LINUX system in PC and SLURM scheduling system
#   on the computing cluster or super-computer.
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
# One need to note that it's just a pseudo-parallism. In other words, this     #
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
#       typing command "chmod +x PACIAE.sh". In addition, one need "make"      #
#       tool. Install "make" using the command, for example on UBUNTU,         #
#       "sudo apt install make".                                               #
#                                                                              #
#   2. Modify the variabls needed.                                             #
#       2.1 Modify the names of source files for Makefile. They shoud be       #
#           the same as the real ones (except for name_x, it's user-defined).  #
#       2.2 Modify the variables: the number of program-running (n_run,        #
#           preferably not larger than the total number of CPU cores of the    #
#           computer), the number of events per program-running (n_eve), and   #
#           other main setups of incident particles and collisions etc.        #
#           "(D=)" means default value.                                        #
#       Here, the variables are aliases for the ones in usu.dat.               #
#       More detailed setups can be modified in the closely following usu.dat  #
#           directly. Search keywords "usu.dat" to find them.                  #
#      2.3 On SLURM/LSF scheduling systems in the computing cluster and the    #
#       sunper-computer. One need modify the additional statements of the      #
#       settings at the beginning of file. (SBATCH, APP_NAME...)               #
#                                                                              #
#   3. Execution.                                                              #
#       3.1 On normal LINUX,  type command "./PACIAE.sh" to run this script.   #
#           The more recommended command is                                    #
#           "time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log", which        #
#                  stores the screen information to a log file.                #
#       3.2 On SLURM: "sbatch PACIAE.sh". On LSF: "bsub PACIAE.sh"             #
#                                                                              #
#                                              By Anke at CCNU on 10/17/2022   #
################################################################################
################################################################################


# Gets current date and time.
echo
current_date_and_time=$(date "+%Y%m%d%H%M%S")
echo "Now is $(date +"%Y-%m-%d %T")"
echo


################################################################################
################################################################################

# Modify the following statements as needed.

################################################################################
# For Makefile
name_x="paciae.x"   # Name of the executable. User-defined.
name_f1="main_23"    # Names of the source files.
name_f2="parini_23"
name_f3="parcas_23"
name_f4="sfm_23"
name_f5="coales_23"
name_f6="hadcas_23"
name_f7="analy"
name_f8="p_23"
# name_f9="eps09"  # Not used now.
# name_f10="xxx"   # Dummy.
################################################################################

################################################################################
# For usu.dat

# Selects the seed of the random number generator in PYTHIA.
iRandom=0   # (D=1) adj(26) selects random generator seed
            #       =0, default PYTHIA seed (19780503), can be used for debug
            #       =1, seed from the real-time system clock (h-min-s-ms)

# Event-number related variables.
n_run=1    # Number of program-running, i.e. total number of CPU cores used.
# n_run=${NP}   # Uncommented this line if submit jods on LSF system.
# ((n_run=SLURM_JOB_NUM_NODES*SLURM_NTASKS_PER_NODE))   # Uncommented for SLURM.
n_eve=10   # neve, number of events per program-running (per core).
n_out=1    # nout, outputs per n_out events.
((tot_eve=n_run*n_eve))   # Total number of events = n_run * n_eve

# Setups of the incident particles.
naProj=197    # nap, nucleon number of projectile
nzProj=79     # nzp, proton number of projectile
naTarg=197    # nat, nucleon number of target
nzTarg=79     # nzt, proton number of target
            # for pp, pA (Ap), and AA
            #  p+p   : 1,1,1,1
            #  p+pbar: 1,1,1,-1
            #  pbar+p: 1,-1,1,1
            #  p+n   : 1,1,1,0
            #  p+p   : 1,0,1,1
            #  n+n   : 1,0,1,0
            #  p+Pb  : 1,1,208,82
            #  Au+Au : 197,79,197,79; Pb+Pb: 208,82,208,82;
            #  Xe+Xe : 129,54,129,54; U + U: 238,92,238,92;
            #  Cu+Cu : 63,29,63,29;   Ag+Ag: 108,47,108,47.
iProj=1     # ipden =0,  if projectile is nucleon (anti-nucleon)
            #       =1,  if nucleus
            #       =2,  for e+e-
            #       =11, if projectile is e- (e+)
            #       =12, if projectile is nu_e (nu_ebar)
            #       =13, if projectile is mu- (mu+)
            #       =14, if projectile is nu_mu (nu_mubar)
            #       =15, if projectile is tau- (tau+)
            #       =16, if projectile is nu_tau (nu_taubar)
iTarg=1     # itden =0,  if target is nucleon (anti-nucleon)
            #       =1,  if nucleus
            #       =2,  for e+e-
            #...
            # for eA, nu_eA, etc.:
            # e^-A:     nap=1, nzp=-1, ipden=11, itden=1,
            # e^+A:     nap=1, nzp=1,  ipden=11, itden=1,
            # nu_eA:    nap=1, nzp=-1, ipden=12, itden=1,
            # nu_ebarA: nap=1, nzp=1,  ipden=12, itden=1,

# Setups of collision.
b_min=8.24     # bmin, min b parameter
b_max=9.23     # bmax, max b parameter
bSamp=2     # (D=2) psno, controls b parameter sampling method.
            #       0=fixed b, 1=systematic, 2=random
iFrame=1    # (D=1) collision frame, 1=collider, 0=fixed target.
cmsE=200    # (GeV) ee, colliding CMS energy / incident momentum

iOverlap=0 # (D=0) adj(30) =0, without more requirements
            #              =1, distributes the participant nucleons in 
            #                  the overlapping region forcely

iChannel=1  # (D=8) nchan, =0: inelastic (INEL)
            #               1: Non Single Difractive (NSD)
            #               2: qqb --> gamma^*/Z^0, used to generate Drell-Yan
            #               3: J/psi production
            #               4: heavy-flavor production
            #               5: direct photon
            #               6: soft only
            #               7: W+/- production
            #               8: pythia default (msel=1)
            #               9: Z0 production
            #               setting of 0,1,3,4,5,7,8 and 9 is ready

iHadcas=1   # (D=1) kjp21, =0, without hadcas; 1, with hadcas
iStage=4    # (D=4) adj(40), =1, stops event after parini
            #                 2, ... parcas
            #                 3, ... 
            #                 4, ... hadcas

# LUND string parameter
aLund=0.3  # (D=0.3)  adj(6), alpha in the LUND string fragmentation function.
bLund=0.58 # (D=0.58) adj(7), beta.

# K factor
kParcas=1           # (D=1) adj(1), K factor multiplying on the differential 
                    #                 cross-section in parton cascade.
                    #                =0, without parton cascade.
kPythia=1           # (D=1.5 / 1) adj(10), K factor multiplying the differential
                    #                      cross sections for hard parton-parton
                    #                      process in PYTHIA.

# String fragmentation and Monte-Carlo coalescence
iHadronization=1    # (D=0) adj(12), model for hadronization
                    #              =0, string fragmentation (sfm)
                    #               1, Monte-Carlo coalescence (coal)
iStringTension=4    # (D=4) kjp22, selects the form of string rension in sfm.
                    #            =1, variable single string tension
                    #             2, variable multiple string tension
                    #             3, variable (single+multiple) string tension
                    #             4, default string tension
iPhase=0    # (D=0) adj(21), 1=with or 0=without phase space constraint in the 
            #                coalescence.
iDeexc=0    # (D=0) adj(29), selects fragmentation function in deexcitation of 
            #       the energetic quark in coalescence
            #       =0, LUND string fragmentation function
            #       =1, Field-Feynman parametrization function
iCME=1     # (D=1) adj(23), choice of chiral magnetic effect(CME)
            #       =0, no CME-induced charge separation mechanism
            #       =1, with CME-induced charge separation mechanism

################################################################################

################################################################################
# The following statements are optional. Usually they are good enough.

# Time accuracy.
dtParini=0.00001    #(D=0.00001) ddt, min distinguishble collision time used in 
                    #                the partonic initiation.
dtParcas=0.03       #(D=0.03) adj(19), time accuracy used in the parton cascade.
dtHadcas=0.1        #(D=0.1) adj(11), time accuracy used in the hadron cascade.

# Statistical histogram intervals in analy. asd(i) in usu.dat.
asd01=0.35   # (D=0.25) y interval
asd02=0.5    # (D=0.5)  pT interval
asd03=0.35   # (D=0.25) eta interval
asd04=0.3    # (D=0.2)  eT interval
asd05=25     # (D=25)   multiplicity interval

pTperturb=0.1  # (D=0.1) smadel, small perturbation of pT ellipse from circle
################################################################################

# More options could be modified on the following "usu.dat" directly.

################################################################################
################################################################################
####################               usu.dat              ########################
echo "${n_eve},${n_out},1                       ! neve,nout,nosc"     > usu.dat
echo "${naProj},${nzProj},${naTarg},${nzTarg}   ! nap,nzp,nat,nzt"   >> usu.dat
echo "${dtParini},0.01,${b_min},${b_max},10     ! ddt,dtt,bmin,bmax,nmax"   >> usu.dat
echo "${iHadcas},${iFrame},1.2,200.,1     ! kjp21,ifram,para7,para10,kjp20" >> usu.dat
echo "3.1416,${iProj},${iTarg}                  ! pio,ipden,itden"   >> usu.dat
echo "20,5,2                                 ! ispmax,isdmax,iflmax" >> usu.dat
echo "211,-211,321,-321,3122,-3122,3312,-3312,3334,-3334   ! KF code of particles" >> usu.dat
echo "2212,-2212,2112,-2112,3212,-3212,3112,3222,310,333   ! KF code of particles" >> usu.dat
echo "${asd01},${asd02},${asd03},${asd04},${asd05}  ! asd(i), i=1,5" >> usu.dat
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
echo "4.5,10.,${cmsE}                           ! parp21,parp22,win" >> usu.dat
echo "0.,.5,1,1,${iChannel}         ! ttaup,taujp,iabsb,iabsm,nchan" >> usu.dat
echo "7.2,4.,${bSamp},40.,20.,0.,0.1  ! para13,para14,psno,para15,para16,ajpsi,vneum" >> usu.dat
echo "40.,40.,25,10.                  ! para1_1,para1_2,para2,para4" >> usu.dat
echo "0.5,0.6,2.,0.05,0.35,20     ! tdh,cptl,cptu,cptl2,cptu2,itnum" >> usu.dat
echo "1,7,0,2,3,2                                   ! mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw" >> usu.dat
echo "${kParcas},0.47,0.4,1000,0,${aLund},${bLund},4,1.9,${kPythia}      ! adj1(1)- adj1(10)  " >> usu.dat
echo "${dtHadcas},${iHadronization},26,18,1.,1,4.,0,${dtParcas},1        ! adj1(11)- adj1(20) " >> usu.dat
echo "${iPhase},4.,${iCME},0.15,0.4,${iRandom},800000.,1.,${iDeexc},${iOverlap}  ! adj1(21)- adj1(30) " >> usu.dat
echo "0.1,0.3,0.4,0.36,1.,0.,100.,3.,2.,${iStage}                        ! adj1(31)- adj1(40) " >> usu.dat
echo "${iStringTension},2,2,0.,0                   ! kjp22,kjp23,kjp24,parp78,mstptj" >> usu.dat
echo "0.,0,${pTperturb},0.05,1.,0.2,0.05           ! parecc,iparres,smadel,dparj4,cp0,cr0,seco" >> usu.dat
echo "0.35776,0.71552                 ! csp_31,csp_32"             >> usu.dat
echo "0.29198,0.58396,0.85464         ! csp_41,csp_42,csp_43"      >> usu.dat
echo "0.23894,0.47788,0.71113,0.91086 ! csp_51,csp_52,csp_53,csp_54" >> usu.dat
echo "0.19964,0.39927,0.59872,0.79707,0.99180 ! csp_61,csp_62,csp_63,csp_64,csp_65" >> usu.dat
####################               usu.dat              ########################
################################################################################
################                                            ####################
################   Annotation is at the end of this file.   ####################
################                                            ####################
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
echo "# src_9 := ${name_f9}"                                        >> Makefile
echo "# scr_10 := ${name_f10}"                                      >> Makefile
echo                                                                >> Makefile
echo "# The intermediate files."                                    >> Makefile
echo "objects := \$(src_1).o \$(src_2).o \$(src_3).o \\"            >> Makefile
echo "           \$(src_4).o \$(src_5).o \$(src_6).o \\"            >> Makefile
echo "           \$(src_7).o \$(src_8).o # \$(src_9).o"             >> Makefile
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
echo "# \$(src_9).o: \$(src_9).f"                                   >> Makefile
echo "# 	\$(compiler) \$(comp_flags) -c -O \$(src_9).f"          >> Makefile
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
    echo
else
    mkdir ./src
fi
# if [ -f "./*.f" ]; then
    mv -f *.f ./src
# fi
mv -f Makefile ./src
wait
# Enters src folder then compiles and builds executable file.
cd ./src 
echo
make
cd ../

# Creates bin folder.
if [ -d "./bin/" ]; then
    echo
else
    mkdir ./bin
fi
cp -f ./src/${name_x} ./bin

# Creates etc folder.
if [ -d "./etc/" ]; then
    echo
else
    mkdir ./etc
fi
# if [ -f "./EPS09*" ]; then
#    mv -f EPS09* ./etc
# fi
mv -f usu.dat ./etc

# Creates log folder.
if [ -d "./log/" ]; then
    echo
else
    mkdir ./log
fi
# Good for normal LINUX and SLURM, but not for LSF.
# TODO(Lei20221016): need to improve for LSF.
mv -f *.out *.log ./log

# File name
if [[ "${naProj}" = "1" && "${nzProj}" = "1" ]]; then
    projName="p"
elif [[ "${naProj}" = "1" && "${nzProj}" = "-1" ]]; then
    projName="pbar"
elif [[ "${naProj}" = "1" && "${nzProj}" = "0" ]]; then
    projName="n"
elif [ "${naProj}" = "63" ]; then
    projName="Cu${naProj}"
elif [ "${naProj}" = "108" ]; then
    projName="Ag${naProj}"
elif [ "${naProj}" = "129" ]; then
    projName="Xe${naProj}"
elif [ "${naProj}" = "197" ]; then
    projName="Au${naProj}"
elif [ "${naProj}" = "208" ]; then
    projName="Pb${naProj}"
elif [ "${naProj}" = "238" ]; then
    projName="U${naProj}"
else
    projName="${naProj}"
fi
# File name
if [[ "${naTarg}" = "1" && "${nzTarg}" = "1" ]]; then
    targName="p"
elif [[ "${naTarg}" = "1" && "${nzTarg}" = "-1" ]]; then
    targName="pbar"
elif [[ "${naTarg}" = "1" && "${nzTarg}" = "0" ]]; then
    targName="n"
elif [ "${naTarg}" = "63" ]; then
    targName="Cu${naTarg}"
elif [ "${naTarg}" = "108" ]; then
    targName="Ag${naTarg}"
elif [ "${naTarg}" = "129" ]; then
    targName="Xe${naTarg}"
elif [ "${naTarg}" = "197" ]; then
    targName="Au${naTarg}"
elif [ "${naTarg}" = "208" ]; then
    targName="Pb${naTarg}"
elif [ "${naTarg}" = "238" ]; then
    targName="U${naTarg}"
else
    targName="${naTarg}"
fi
# File name
if [ "${iFrame}" = "1" ]; then
    collName="COLL"
elif [ "${iFrame}" = "0" ]; then
    collName="FXT"
else
    collName="${iFrame}"
fi
# Creates files.
if [ -d "./${projName}_${targName}_${collName}/" ]; then
    echo
else
    mkdir ./${projName}_${targName}_${collName}
fi
cd ./${projName}_${targName}_${collName}

if [ -d "./${cmsE}GeV/" ]; then
    echo
else
    mkdir ./${cmsE}GeV
fi
cd ./${cmsE}GeV
#ls -a

# b: impact parameter; 
# aL: aLund
# bL: bLund
# K : kPythia (in parini)
#     kI: kPythia (in parini)
#     kC: kParcas
# iO: iOverlap
# iH: iHadronization
# iHc: iHadcas
# iStg: iStage
# iC : iChannel
# iStr: iStringTension
# iP : iPhase
# iD : iDeexc
# iM : icMe
# Default 
dir="./b${b_min}_${b_max}_aL${aLund}_bL${bLund}_K${kPythia}_iStg${iStage}"
# dir="./b${b_min}_${b_max}_aL${aLund}_bL${bLund}_kI${kPythia}_kC${kParcas}_iO${iOverlap}_iH${iHadronization}_iHc${iHadcas}_iStr${iStringTension}_iP${iPhase}_iD${iDeexc}_iM${iCME}_iC${iChannel}_iStg${iStage}"

if [[ "${naProj}" = "1" && "${naTarg}" = "1" ]]; then
# Elementary NN collision
    if [[ "${iHadronization}" = "0" ]]; then
        dir="./Sfm_aL${aLund}_bL${bLund}_K${kPythia}_iStg${iStage}"
    elif [[ "${iHadronization}" = "1" ]]; then
        dir="./Coal_K${kPythia}_iStg${iStage}"
    # else
    fi
# elif [[  ]]; then
else
# pA, AA collisions
    if [[ "${iHadronization}" = "0" ]]; then
        dir="./b${b_min}_${b_max}_sfm_aL${aLund}_bL${bLund}_K${kPythia}_iStg${iStage}"
    elif [[ "${iHadronization}" = "1" ]]; then
        dir="./b${b_min}_${b_max}_coal_K${kPythia}_iStg${iStage}"
    # else
    fi
fi


if [ -d "${dir}/" ]; then
    echo
else
    mkdir ${dir}
fi
#ls -a

if [ -d "${dir}/${tot_eve}_events/" ]; then
    echo
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
# neve,nout,nosc
# neve: events number to be generated
# nout: output the event per nout events
# nosc: OSCAR stander output (oscar.out)
#
# nap,nzp,nat,nzt
# for pA, AA, etc.
# nap(nzp): nucleons (protons) number of projectile
# nat(nzt): nucleons (protons) number of target
#
# for eA, nu_eA, etc.
# e^-A: nap=1,nzp=-1,ipden=11,itden=1, 
# e^+A: nap=1,nzp=1,ipden=11,itden=1, 
# nu_eA: nap=1,nzp=-1,ipden=12,itden=1,
# nu_ebarA: nap=1,nzp=1,ipden=12,itden=1.
#
# ddt,dtt,bmin,bmax,nmax
#  ddt: minimum distinguishble collision time interval used in partonic 
#       initiation in parini.f
#  dtt: time accuracy (not used now)
#  bmin: minimum impact parameters, 
#  bmax: maximum impact parameters,
#  nmax: the number of intervals segmented in [bmin,bmax] when psno=1
#
# kjp21,ifram,para7,para10,kjp20
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
# pio,ipden,itden 
#  pio: pi=3.1416
#  ipden: =0, if projectile is proton (anti-proton)
#         =2, for e+e-
#         =11, if projectile is e- (e+)
#         =12, if projectile is nu_e (nu_ebar)
#         =13, if projectile is mu- (mu+)
#         =14, if projectile is nu_mu (nu_mubar)
#         =15, if projectile is tau- (tau+)
#         =16, if projectile is nu_tau (nu_taubar)
#  itden: =0, if target is proton (anti-proton)
#         =2, for e+e-
#
# ispmax,isdmax,iflmax
#  ispmax: maximum # of different particle pieces to be considered    
#  isdmax: maximum # of different distributions to be considered
#  iflmax: maximum # of windows to be set, =0 means no window at all
#
# ispkf(i,i=1,ispmax):
# KF code: particle code used in PYTHIA and PACIAE, 
#          (see detail in reference: arXiv:hep-ph/0603175) 
#
# asd(i=1,isdmax): interval of the i-th distribution
#  for pp, pbarp etc.
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
# afl(j,i,2): upper-boundary of i-th window for j-th particle
#  for pp, pbarp etc.
#      i=1, rapidity window
#       =2, transverse monmentum 
#       =3, pesudorapidity 
#  for ep, nu_ep, etc. 
#      i=1, Q^2=-q^2 window
#       =2, W^2 
#       =3, y 
#       =4, p_h (haron momentum)
#       =5: z 
#
# parp21,parp22,ee  
#  parp21: lowest CM energy running 'pythia' if nchan=6
#  parp22: lowest CM energy for running 'pythia' if nchan=3
#  ee= cms energy if ifram=1 (collider)
#    = incident momentum if ifram=0 (fixed target)
#
# ttaup,taujp,iabsb,iabsm,nchan
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
# para13,para14,psno,para15,para16,ajpsi,vneum
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
# para1_1,para1_2,para2,para4  
#  para1_1: nn total cross section used in parton initiation
#  para1_2: nn total cross section used in hadron cascade
#  para2: totle cross-section of pi-nucleon 
#  para4: totle cross-section of pi-pi
#
# tdh,cptl,cptu,cptl2,cptu2,itnum
#  tdh: time step used in subroutine 'flow'
#  cptl: lower pt cut in 'flow' for particle 1
#  cptu: upper pt cut in 'flow' for particle 1
#  cptl2: lower pt cut in 'flow' for particle 2
#  cptu2: upper pt cut in 'flow' for particle 2
#  itnum: number of time steps used in 'flow'
#
# mstu21,mstj1_1,mstj1_2,mstj2,mstj3,itorw
#  mstu21: parameter mstu(21) in PYTHIA
#  mstj1_1: =0, no jet fragmentation set in "parini.f"
#  mstj1_2: =1, Lund string fragmentation set in "sfm.f"
#  mstj1_2: =2, independent fragmentation set in "sfm.f"
#  mstj2: parameter mstj(2) in PYTHIA
#  mstj3: parameter mstj(3) in PYTHIA
#  itorw: =1 executing pyevnt, =2 executing pyevnw
#
# adj1(i), i=1,40: switches and/or parameters
#     i=1: K factor in parton cascade
#       2: alpha_s, effective coupling constant in parton rescattering
#       3: parameter \mu^2 (tcut in program), the regulation factor introduced 
#          in the parton-parton differential cross section
#       4: parameter `idw', the number of intervals in the numerical integration
#       5: =0, w/o nuclear shadowing,
#          =1, Wang's nuclear shadowing (PLB 527(2002)85).
#       6: alpha in the LUND string fragmentation function (parj(41) in PYTHIA)
#       7: beta in the LUND string fragmentation function (parj(42) in PYTHIA)
#       8: mstp(82) in PYTHIA 
#       9: parp(81) in PYTHIA
#       10: K factor (parp(31) in PYTHIA).
#       11: time accuracy used in the hadron cascade
#       12: model for hadronization 
#           =0 string fragmentation,
#           =1 Monte Carlo coalescence model. 
#       13: dimension of meson table considered in coalescence model.
#       14: dimension of baryon table considered coalescence model.
#       15: string tension of qqbar simple string
#       16: number of loops in the deexcitation of energetic quark in the 
#           Monte Carlo coalescence model.
#       17: the threshold energy in the deexcitation of energetic quark in 
#           the Monte Carlo coalescence model.
#       18: with or without Pauli blocking in the parton cascade
#           =0, without Pauli blocking
#           =1, with Pauli blocking.
#       19: time accuracy used in the parton cascade 
#       20: the optional parton-parton cross section in the parton rescattering
#           =0, LO pQCD parton-parton cross section 
#           =1, keeping only leading divergent terms in the LO pQCD 
#              parton-parton cross section (B. Zhang)
#           =2, the same as 0 but flat scattering angle distribution is assumed
#           =3, the same as 1 but flat scattering angle distribution is assumed.
#       21: with or without phase space constraint in the Monte Carlo 
#           coalescence model
#          =0, without phase space constraint
#          =1, with phase space constraint
#       22: critical value (D=$) of the product of radii in position and 
#           momentum phase spaces 
#       23: choice of chiral magnetic effect(CME)
#           =0: no CME-induced charge separation mechanism
#           =1: with CME-induced charge separation mechanism
#       24: the virtuality cut ('tl0' in program) in the time-like radiation in 
#           parton rescattering
#       25: $\Lambda_{QCD}$ in PYTHIA
#       26: selection of random number seed,
#           =0, default PYTHIA seed (19780503), can be used for debug
#           =other, seed from the real-time clock
#       27: largest momentum allowed for produced particle
#       28: concerned to the largest position allowed for produced particle
#       29: =0, Lund string fragmentation function is used in coalescence model
#           =1, Field-Feymman fragmentation function is used in coalescence model
#       30: =0, without more requirements
#           =1, distributes the participant nucleons in overlapping areas forcely
#       31: parj(1) in PYTHIA 
#       32: parj(2) in PYTHIA 
#       33: parj(3) in PYTHIA 
#       34: parj(21) in PYTHIA 
#       35: mstp(91) in PYTHIA, parton transverse momentum ($k_{\perp}$) 
#           distribution inside hadron
#           =1: Gaussian
#           =2: exponential
#       36: with or without phenomenological parton energy loss in parton 
#           rescattering
#           =0, without 
#           =1, with 
#       37: the coefficient in phenomenological parton energy loss
#       38: $p_T$ cut in phenomenological parton energy loss
#       39: width of Gaussian parton k_perp distribution in hadron if mstp(91)=1
#           width of exponential k_{\perp} distribution in hadron if mstp(91)=2
#       40: optional event stopping point 
#          =1, after parton initiation,
#          =2, after parton rescattering,
#          =4, after hadron rescattering 
#
# kjp22,kjp23,kjp24,parp78,mstptj
#  kjp22: =1, variable single string tension and parj(1) etc.
#         =2, variable multiple string tension and parj(1) etc.
#         =3, variable (single+multiple) string tension and parj(1) etc.
#         =4, default string tension and parj(1) etc.
#         see PhysRevC.98.034917(2018).
#  kjp23: optional model for the calculation of participant nucleon # (npart)
#         =1, geometric model
#         =2, Glauber model
#  kjp24: optional distribution in Glauber model
#         =1, sharp sphere
#         =2, Woods-Saxon
#  parp78: parameter controling amount of colour reconnection in final state
#  mstptj: input value of mstp(111) (mstj(1)) for pp, pA (AP), and AA 
#         (for e+e-), i.e. switch off/on (0/1) hadronization in PYTHIA (parini).
#
# parecc,iparres,smadel,dparj4,cp0,cr0,seco
#  parecc: proportional factor between initial spatial space eccentricity and 
#          final momentum space ellipticity 
#  iparres: =0 consider elastic parton-parton collisions only in parton 
#              rescattering
#           =1 otherwise
#  smadel: small perturbation of ellipse from circle
#  dparj4: default parj(4)
#  cp0,cr0: parameters in parameterization of multiple string effect
#  seco: parameter in popcorn mechanism for correction of parj(1)
#
# cumulent sum of flavor generation probability, for an example: 
# p_31=(1-amd/sms)/2., p_32=(1-amu/sms)/2., p_33=(1-ams/sms)/2.; 
# amd (amu, ams): d (u,s) quark constituent mass, sms=amd+amu+ams; 
# p_31+p_32+p33=1;
# cumulent sum: csp_31=p_31, csp_32=p_31+p_32
# u and d have same constituent mass and probability
######################        Annotation of usu.dat         ####################
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
