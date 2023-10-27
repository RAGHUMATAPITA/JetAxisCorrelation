#!/bin/bash

#root -l -b -q 'Tree_Analyzer.C("inputfile_PbPb.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/Hydjet_pThat15", "PbPb", 1)'
root -l -b -q 'Tree_Analyzer.C("inputfile_PbPb_MAODData.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/Data/rootfiles", "PbPb", 0)'


#root -l -b -q 'Tree_Analyzer.C("inputfile_pp.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/Pythia_pThat15", "pp", 1)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_pprefData.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/Data/rootfiles", "pp", 0)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_pprefData_New.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/Data/rootfiles", "pp", 0)'
