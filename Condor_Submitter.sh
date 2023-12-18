#!/bin/bash

rm -rf *.C~
rm -rf *.h~
rm -rf *.sh~
rm -rf *.txt~

running_system=$4

file="Data_MC_txtFiles/Pythia_pThat15.txt"

if [ $running_system == 1 ]; then
    #file="Data_MC_txtFiles/PbPb_2018_AODMC_New_TIFR_Path.txt"
    file="Data_MC_txtFiles/Hydjet_pThat15.txt"
    echo ${file}
elif [ $running_system == 2 ]; then
    #file="Data_MC_txtFiles/ppref_2017_AODMC_New_TIFR_Path.txt"
    file="Data_MC_txtFiles/Pythia_pThat15.txt"
    echo ${file}
elif [ $running_system == 3 ]; then
    file="Data_MC_txtFiles/PbPb_2018_AODData_New_TIFR_Path.txt"
    echo ${file}
else 
    file="Data_MC_txtFiles/ppref_2017_AODData_New_TIFR_Path.txt"
    echo ${file}
fi

istestrun=$5
testRun_condition=true

if [ $istestrun == 1 ]; then
    testRun_condition=true
else
    testRun_condition=false
fi

isMC=$6
mc_or_data=1

if [ $isMC == 1 ]; then
    mc_or_data=1 
else 
    mc_or_data=0
fi
#job flavour (espresso=20min, microcentury=1h, longlunch=2h, workday=8h, tomorrow=1d, testmatch=3d, nextweek=1w)
# Arguments: 1 --> input txt file, 2 --> job flavour, 3 --> number of cpus, 4 --> outfile txt file, 5 --> path to save root files

if [ $istestrun == 1 ]; then
    if [ $isMC == 1 ]; then
	root -l -b -q "file_write_copy.C(${testRun_condition}, \"${file}\", \"$1\", \"espresso\", 4, \"/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles\", \"$2\", \"$3\", ${mc_or_data})"
    else
	root -l -b -q "file_write_copy.C(${testRun_condition}, \"${file}\", \"$1\", \"espresso\", 4, \"/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/Data/rootfiles\", \"$2\", \"$3\", ${mc_or_data})"
    fi
else
    if [ $isMC == 1 ]; then
	root -l -b -q "file_write_copy.C(${testRun_condition}, \"${file}\", \"$1\", \"workday\", 4, \"/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles\", \"$2\", \"$3\", ${mc_or_data})"
    else 
	root -l -b -q "file_write_copy.C(${testRun_condition}, \"${file}\", \"$1\", \"workday\", 4, \"/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/Data/rootfiles\", \"$2\", \"$3\", ${mc_or_data})"
    fi
fi

condor_submit Submission_Files/$1/$2/Submission_$1_$2.sub
