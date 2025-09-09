# -*- coding: utf-8 -*-
"""
This code is part of the Pipeline-Codes
Part 2. Data Processing Part 1. Data Normalization, Reduction

Input: 
    Type : XX_Balanced_Training.txt

Read carefully the ParamFile (Step2_Parameters.txt)    
    
Processing Format:
    python Step1SmilesBalancingPartition.py -i InputFile -p ParamFile
    python Step1SmilesBalancingPartition.py -h
"""

from Modelling_Lib import  SimplestDataReduction
import argparse
import sys
import numpy as np
from sklearn import preprocessing
import pickle


def ReadingParamFile (Filename):
    f = open(Filename, "r")
    FinalListElements = 0
    for line in f:
        if line != "\n" and line[0] != "#":
            Temp = line.split("\n")[0].split(" = ")
            if Temp[0] == "Data_Normalization":
                if Temp[1]=="True":
                    Data_Normalization=True
                    FinalListElements=FinalListElements+1
                else:
                    if Temp[1]=="False":
                        Data_Normalization=False
                        FinalListElements=FinalListElements+1
                    else:
                        print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                        return None
            if Temp[0] == "Simple_Variable_Reduction":
                if Temp[1]=="True":
                    Simple_Variable_Reduction=True
                    FinalListElements=FinalListElements+1
                else:
                    if Temp[1]=="False":
                        Simple_Variable_Reduction=False
                        FinalListElements=FinalListElements+1
                    else:
                        print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                        return None

            if Temp[0] == "Var_Curoff":
                try:
                    Var_Curoff=float(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None

    f.close()
    return Data_Normalization, Simple_Variable_Reduction, Var_Curoff




parser = argparse.ArgumentParser(description="Step 2 of the Pipeline: Data Curation Part 1.")
parser.add_argument("-i", type=str, help="Input file XXX_Balanced_Training.txt from Step 1")
parser.add_argument("-p", type=str, help="Input file with computation parameters.")


if len(sys.argv) == 1:  # No arguments provided
    sys.argv = [
        'Step2DataProcessing1 1.py',
        '-i', 'HCT-15_Curated_Balanced_Training.txt',  # Replace with your actual input file
        '-p', 'Step2ParamFile.txt'   # Replace with your parameter file
    ]

args = parser.parse_args()

if args.i == None:
    print("Missing Input file")
    sys.exit()
if args.p == None:
    print("Missing Parameter file")
    sys.exit()

#Reading Parameter File
ParamFile = args.p
print(ReadingParamFile (ParamFile))
Data_Normalization, Simple_Variable_Reduction, Var_Curoff = ReadingParamFile (ParamFile)

if Data_Normalization==False:
    if Simple_Variable_Reduction==False:
        print("Both, normalization and Variable Reduction are set to False, the Step 2 is nor requiered. Exit")
        sys.exit()
        

#Reading Input File
InputFile = args.i
f = open(InputFile, "r")
glog = open(InputFile.split(".")[0]+"_Step2_logFile.log", "w")
X=[]
N = 0
for line in f:
    if N > 0:
        Temp = line.split("\n")[0].split("\t")
        Xtemp = Temp[3:]
        TypeTemp = " ".join(Xtemp)
        if TypeTemp.find(".")>-1:
            Xtemp = [float(a) for a in Xtemp]
        else:
            Xtemp = [int(a) for a in Xtemp]
        X.append(Xtemp)
    else:
        FPAll = line.split("\n")[0].split("\t")
        FPAll = FPAll[3:]
    N=N+1
f.close()

X = np.array(X)

if Data_Normalization:
    glog.write("############################################################################\n")
    glog.write("#### BEGIN DATA NORMALIZATION / VARIABLE SELECTION #######\n")
    scaler = preprocessing.StandardScaler().fit(X)
    Step2Saved=[scaler]
    glog.write("# The normilized scaler is located at " +InputFile.split(".")[0]+"_Step2_Scaler.pkl"+ "\n")
    if Simple_Variable_Reduction:
        Selected = SimplestDataReduction(X, Var_Curoff)
        FPselected = [FPAll[i] for i in Selected]
        Step2Saved = Step2Saved + [FPselected]
        glog.write("# The selected variables are located at " +InputFile.split(".")[0]+"_Step2_Scaler.pkl"+ "\n")
        glog.write("# The selected variables are also presented here:"+ "\n")
        for ss in FPselected:
            glog.write(ss + "\n")
    else:
        Step2Saved = Step2Saved + [None]
        glog.write("# No variable section was included in the Step 2 analysis \n")
else:
    glog.write("############################################################################\n")
    glog.write("#### BEGIN DATA NORMALIZATION / VARIABLE SELECTION #######\n")
    glog.write("# No variable normalization section was included in the Step 2 analysis \n")
    if Simple_Variable_Reduction:
        Selected = SimplestDataReduction(X, Var_Curoff)
        FPselected = [FPAll[i] for i in Selected]
        Step2Saved = [None, FPselected]
        glog.write("# The selected variables are located at " +InputFile.split(".")[0]+"_Step2_Scaler.pkl"+ "\n")
        glog.write("# The selected variables are also presented here:"+ "\n")
        for ss in FPselected:
            glog.write(ss + "\n")
    else:
        Step2Saved = [None, None]
        glog.write("# No variable section was included in the Step 2 analysis \n")

with open(InputFile.split(".")[0]+"_Step2_Scaler.pkl",'wb') as f:
    pickle.dump(Step2Saved,f)
        
        
glog.write("#### END DATA NORMALIZATION / VARIABLE SELECTION #######\n")
glog.write("############################################################################\n")
glog.write("SUCCESS COMPLETED\n")

glog.close()








