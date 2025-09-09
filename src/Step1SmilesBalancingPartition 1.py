# -*- coding: utf-8 -*-
"""
This code is part of the Pipeline-Codes
Part 1. Smiles Processing and Balancing

Input: 
    Type 1: ID/Smiles/Class, with Heather

Read carefully the ParamFile (Step1_Parameters.txt)    
    
Processing Format:
    python Step1SmilesBalancingPartition.py -i InputFile -p ParamFile
    python Step1SmilesBalancingPartition.py -h
"""


from Modelling_Lib import  SmilesSanitaze, ComputeDescriptors, DataBalancing
import argparse
import sys
import numpy as np
from sklearn.model_selection import train_test_split


def ReadingParamFile (Filename):
    f = open(Filename, "r")
    FinalListElements = 0
    for line in f:
        if line != "\n" and line[0] != "#":
            Temp = line.split("\n")[0].split(" = ")
            if Temp[0] == "Fix_Salts":
                if Temp[1]=="True":
                    IgnoreSalts=True
                    FinalListElements=FinalListElements+1
                else:
                    if Temp[1]=="False":
                        IgnoreSalts=False
                        FinalListElements=FinalListElements+1
                    else:
                        print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                        return None
            if Temp[0] == "Include_Isomers":
                if Temp[1]=="True":
                    IncludeIsomer=True
                    FinalListElements=FinalListElements+1
                else:
                    if Temp[1]=="False":
                        IncludeIsomer=False
                        FinalListElements=FinalListElements+1
                    else:
                        print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                        return None
            if Temp[0] == "Ignore_Charge":
                if Temp[1]=="True":
                    IgnoreCharge=True
                    FinalListElements=FinalListElements+1
                else:
                    if Temp[1]=="False":
                        IgnoreCharge=False
                        FinalListElements=FinalListElements+1
                    else:
                        print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                        return None
            if Temp[0] == "Descriptors_Type":
                try:
                    FingerPrintType=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "Partition_Ratio":
                try:
                    PartitionRatio=float(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "Metric_Cluster":
                try:
                    MetricCluster=Temp[1]
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "N_Cluster_Min":
                try:
                    Cmin=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "Split_Ratio":
                try:
                    SplitRatio=Temp[1].split("-")
                    SplitRatio=[int(a) for a in SplitRatio]
                    if sum(SplitRatio)==100:
                        FinalListElements=FinalListElements+1
                    else:
                        print("Check the sum in " + Temp[0])
                        return None
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
    f.close()
    return IgnoreSalts, IncludeIsomer, IgnoreCharge, FingerPrintType, PartitionRatio, MetricCluster, Cmin, SplitRatio




parser = argparse.ArgumentParser(description="Step 1 of the Pipeline: Smiles curation and data processing.")
parser.add_argument("-i", type=str, help="Input file with the structure ID \t Smile \t Class. First line is the header")
parser.add_argument("-p", type=str, help="Input file with computation parameters.")
args = parser.parse_args()

#if args.i == None:
 #   print("Missing Input file")
  #  sys.exit()
#if args.p == None:
 #   print("Missing Parameter file")
  #  sys.exit()

#Reading Parameter File
ParamFile = "Step1ParamFile.txt"
print(ReadingParamFile (ParamFile))
IgnoreSalts, IncludeIsomer, IgnoreCharge, FingerPrintType, PartitionRatio, MetricCluster, Cmin, SplitRatio = ReadingParamFile (ParamFile)

#Reading Input File
InputFile = "HT-29.txt"
f = open(InputFile, "r")
g = open(InputFile.split(".")[0]+"_Curated.txt", "w")
glog = open(InputFile.split(".")[0]+"_Step1_logFile.log", "w")
SmileID, SmilesF, Class, SmilesDescription = [],[],[], []
LinesError = []
N = 0
Nerror =0

"""
print("Starting SMILES processing...")

for line in f:
    if N > 0:
        # Progress monitoring (console only)
        Temp = line.split("\n")[0].split("\t")
        
        if N % 100 == 0:
            print(f"Processing compound {N}/13016 ({N/130.16:.1f}%)")
        
        if N % 2000 == 0:
            print(f"Current compound: {Temp[0]}")
        
        
        SmileCorrected = SmilesSanitaze(Temp[1], IgnoreSalts, IncludeIsomer, IgnoreCharge)
        if SmileCorrected != None:
            DescriptorsT= ComputeDescriptors (SmileCorrected, FingerPrintType, IncludeIsomer)
            if DescriptorsT != None:
                g.write(Temp[0] + "\t" + SmileCorrected + "\t" + Temp[2]+ "\n")
                SmileID.append(Temp[0])
                SmilesF.append(SmileCorrected)
                Class.append(int(Temp[2]))
                SmilesDescription.append(DescriptorsT)
            else:
                Nerror = Nerror+1
                LinesError.append(line)
        else:
            Nerror = Nerror+1
            LinesError.append(line)
    else:
        g.write(line)
    N=N+1
f.close()
g.close()

print("SMILES processing completed!")

# EXACT LOG FORMAT - BEGIN EXCLUDED SMILES
if Nerror>0:
    glog.write("############################################################################\n")
    glog.write("#### BEGIN EXCLUDED SMILES #######\n")
    glog.write("A total of " + str(N-1) + " smiles were processed" + "\n")
    glog.write("The full curated table is presented in the file: "+  InputFile.split(".")[0]+"_Curated.txt"+ "\n")
    glog.write("A total of " + str(Nerror) + " smiles were excluded and showed bellow" + "\n")
    for line in LinesError:
        glog.write(line)
    glog.write("#### END EXCLUDED SMILES #######\n")
    glog.write("############################################################################\n")
else:
    glog.write("############################################################################\n")
    glog.write("#### BEGIN EXCLUDED SMILES #######\n")
    glog.write("A total of " + str(N-1) + " smiles were processed" + "\n")
    glog.write("The full curated table is presented in the file: "+  InputFile.split(".")[0]+"_Curated.txt"+ "\n")
    glog.write("No smiles were excluded" + "\n")
    glog.write("#### END EXCLUDED SMILES #######\n")
    glog.write("############################################################################\n")

# EXACT LOG FORMAT - BEGIN BALANCING DATA
glog.write("############################################################################\n")
glog.write("#### BEGIN BALANCING DATA #######\n")

Np = np.sum(Class)
Nn = len(Class)-Np
ToBalance = False

if Np>=Nn:
    if Nn/Np < PartitionRatio:
        ToBalance = True
        FractionComputed = Nn/Np
        MajorClass = "1"
else:
    if Np/Nn < PartitionRatio:
        ToBalance = True
        FractionComputed = Np/Nn
        MajorClass = "0"
        
if ToBalance:
    print("🔄 Data balancing required - finding optimal clusters...")
    print("This may take 1-4 hours for large datasets...")
    
    SmiID = []
    SelectedCases, NC = DataBalancing (SmilesDescription, Class, MetricCluster, Cmin)
    
    print(f"✅ Clustering completed! Optimal clusters found: {NC}")
    
    # EXACT LOG FORMAT - DATA NEED TO BE BALANCED
    glog.write("#### DATA NEED TO BE BALANCED #######\n")
    glog.write("Computed fraction is " + str(FractionComputed) + "\n")
    glog.write("Majoritary Class is " + MajorClass + "\n")
    glog.write("Number of optimal clusters is " + str(NC) + "\n")
    glog.write("The full balanced table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced.txt"+ "\n")
    glog.write("The exluded clases after balancing are in the file: "+  InputFile.split(".")[0]+"_Curated_Excluded.txt"+ "\n")
    glog.write("#### END BALANCING DATA #######\n")
    glog.write("############################################################################\n")
    
    FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
    FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
    f = open(InputFile.split(".")[0]+"_Curated_Balanced.txt", "w")
    f.write("\t".join(FP) + "\n")
    g = open(InputFile.split(".")[0]+"_Curated_Excluded.txt", "w")
    g.write("SMILEID"+"\n")
    k=0
    for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
        if k in SelectedCases:
            FP = [str(k) for k in smifp]
            FP = "\t".join(FP)
            f.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
            SmiID.append(sid)
        else:
            g.write(sid + "\n")
        k=k+1
    f.close()
    g.close()
else:
    print("✅ Data is already balanced - no clustering needed")
    SmiID = []
    # EXACT LOG FORMAT - DATA DONT NEED TO BE BALANCED
    glog.write("#### DATA DONT NEED TO BE BALANCED #######\n")
    glog.write("Computed fraction is " + str(FractionComputed) + "\n")
    glog.write("Majoritary Class is " + MajorClass + "\n")
    glog.write("The full table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced.txt"+ "\n")
    glog.write("#### END BALANCING DATA #######\n")
    glog.write("############################################################################\n")
    
    FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
    FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
    f = open(InputFile.split(".")[0]+"_Curated_Balanced.txt", "w")
    f.write("\t".join(FP) + "\n")
    k=0
    for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
        FP = [str(k) for k in smifp]
        FP = "\t".join(FP)
        f.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        SmiID.append(sid)
    f.close()

print("🔄 Starting data partitioning...")

# EXACT LOG FORMAT - BEGIN DATA PARTITION
glog.write("############################################################################\n")
glog.write("#### BEGIN DATA PARTITION #######\n")

SmiID = np.array(SmiID)

if len(SplitRatio)==3:    
    glog.write("The Training table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_Training.txt"+ "\n")
    glog.write("The Test table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_Test.txt"+ "\n")
    glog.write("The External table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_External.txt"+ "\n")
    
    NTest = int((SplitRatio[1]/100)*len(SelectedCases))
    NExt = int((SplitRatio[2]/100)*len(SelectedCases))
    NTrain = len(SelectedCases) - NTest - NExt
    SmiIDt, SmiID_ext = train_test_split(SmiID, test_size = NExt, random_state=42)
    SmiID_train, SmiID_test = train_test_split(SmiIDt, test_size = NTest, random_state=42)
    
    FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
    FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
    f1 = open(InputFile.split(".")[0]+"_Curated_Balanced_Training.txt", "w")
    f1.write("\t".join(FP) + "\n")
    f2 = open(InputFile.split(".")[0]+"_Curated_Balanced_Test.txt", "w")
    f2.write("\t".join(FP) + "\n")
    f3 = open(InputFile.split(".")[0]+"_Curated_Balanced_External.txt", "w")
    f3.write("\t".join(FP) + "\n")

    for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
        FP = [str(k) for k in smifp]
        FP = "\t".join(FP)
        if sid in SmiID_train:
            f1.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        if sid in SmiID_test:
            f2.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        if sid in SmiID_ext:
            f3.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
    f1.close()
    f2.close()
    f3.close()
else:
    if len(SplitRatio)==2:    
        glog.write("The Training table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_Training.txt"+ "\n")
        glog.write("The External table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_External.txt"+ "\n")
        NExt = int((SplitRatio[1]/100)*len(SelectedCases))
        SmiID_train, SmiID_ext = train_test_split(SmiID, test_size = NExt, random_state=42)
        FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
        FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
        f1 = open(InputFile.split(".")[0]+"_Curated_Balanced_Training.txt", "w")
        f1.write("\t".join(FP) + "\n")
        f2 = open(InputFile.split(".")[0]+"_Curated_Balanced_External.txt", "w")
        f2.write("\t".join(FP) + "\n")

        for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
            FP = [str(k) for k in smifp]
            FP = "\t".join(FP)
            if sid in SmiID_train:
                f1.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
            if sid in SmiID_ext:
                f2.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        f1.close()
        f2.close()
        
    else:
        print("An error ocurred during data partition, please check the Param File")
        sys.exit()

# EXACT LOG FORMAT - END DATA PARTITION
glog.write("#### END DATA PARTITION #######\n")
glog.write("############################################################################\n")
glog.write("SUCCESS COMPLETED\n")

print("🎉 SUCCESS COMPLETED!")
print("Check the log file for detailed results.")

glog.close()

"""

for line in f:
    if N > 0:
        Temp = line.split("\n")[0].split("\t")
        SmileCorrected = SmilesSanitaze(Temp[1], IgnoreSalts, IncludeIsomer, IgnoreCharge)
        if SmileCorrected != None:
            DescriptorsT= ComputeDescriptors (SmileCorrected, FingerPrintType, IncludeIsomer)
            if DescriptorsT != None:
                g.write(Temp[0] + "\t" + SmileCorrected + "\t" + Temp[2]+ "\n")
                SmileID.append(Temp[0])
                SmilesF.append(SmileCorrected)
                Class.append(int(Temp[2]))
                SmilesDescription.append(DescriptorsT)
            else:
                Nerror = Nerror+1
                LinesError.append(line)
        else:
            Nerror = Nerror+1
            LinesError.append(line)
    else:
        g.write(line)
    N=N+1
f.close()
g.close()

if Nerror>0:
    glog.write("############################################################################\n")
    glog.write("#### BEGIN EXCLUDED SMILES #######\n")
    glog.write("A total of " + str(N) + " smiles were processed" + "\n")
    glog.write("The full curated table is presented in the file: "+  InputFile.split(".")[0]+"_Curated.txt"+ "\n")
    glog.write("A total of " + str(Nerror) + " smiles were excluded and showed bellow" + "\n")
    for line in LinesError:
        glog.write(line)
    glog.write("#### END EXCLUDED SMILES #######\n")
    glog.write("############################################################################\n")
else:
    glog.write("############################################################################\n")
    glog.write("#### BEGIN EXCLUDED SMILES #######\n")
    glog.write("A total of " + str(N) + " smiles were processed" + "\n")
    glog.write("The full curated table is presented in the file: "+  InputFile.split(".")[0]+"_Curated.txt"+ "\n")
    glog.write("No smiles were excluded" + "\n")
    glog.write("#### END EXCLUDED SMILES #######\n")
    glog.write("############################################################################\n")



#Balancing
glog.write("############################################################################\n")
glog.write("#### BEGIN BALANCING DATA #######\n")

Np = np.sum(Class)
Nn = len(Class)-Np
ToBalance = False

if Np>=Nn:
    if Nn/Np < PartitionRatio:
        ToBalance = True
        FractionComputed = Nn/Np
        MajorClass = "1"
else:
    if Np/Nn < PartitionRatio:
        ToBalance = True
        FractionComputed = Np/Nn
        MajorClass = "0"
        
if ToBalance:
    SmiID = []
    SelectedCases, NC = DataBalancing (SmilesDescription, Class, MetricCluster, Cmin)
    glog.write("#### DATA NEED TO BE BALANCED #######\n")
    glog.write("Computed fraction is " + str(FractionComputed) + "\n")
    glog.write("Majoritary Class is " + MajorClass + "\n")
    glog.write("Number of optimal clusters is " + str(NC) + "\n")
    glog.write("The full balanced table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced.txt"+ "\n")
    glog.write("The exluded clases after balancing are in the file: "+  InputFile.split(".")[0]+"_Curated_Excluded.txt"+ "\n")
    glog.write("#### END BALANCING DATA #######\n")
    glog.write("############################################################################\n")
    
    FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
    FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
    f = open(InputFile.split(".")[0]+"_Curated_Balanced.txt", "w")
    f.write("\t".join(FP) + "\n")
    g = open(InputFile.split(".")[0]+"_Curated_Excluded.txt", "w")
    g.write("SMILEID"+"\n")
    k=0
    for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
        if k in SelectedCases:
            FP = [str(k) for k in smifp]
            FP = "\t".join(FP)
            f.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
            SmiID.append(sid)
        else:
            g.write(sid + "\n")
        k=k+1
else:
    SmiID = []
    glog.write("#### DATA DONT NEED TO BE BALANCED #######\n")
    glog.write("Computed fraction is " + str(FractionComputed) + "\n")
    glog.write("Majoritary Class is " + MajorClass + "\n")
    glog.write("The full table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced.txt"+ "\n")
    glog.write("#### END BALANCING DATA #######\n")
    glog.write("############################################################################\n")
    FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
    FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
    f = open(InputFile.split(".")[0]+"_Curated_Balanced.txt", "w")
    f.write("\t".join(FP) + "\n")
    k=0
    for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
        FP = [str(k) for k in smifp]
        FP = "\t".join(FP)
        f.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        SmiID.append(sid)



#Data partition and normalization

glog.write("############################################################################\n")
glog.write("#### BEGIN DATA PARTITION #######\n")

SmiID = np.array(SmiID)

if len(SplitRatio)==3:    
    glog.write("The Training table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_Training.txt"+ "\n")
    glog.write("The Test table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_Test.txt"+ "\n")
    glog.write("The External table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_External.txt"+ "\n")
    
    NTest = int((SplitRatio[1]/100)*len(SelectedCases))
    NExt = int((SplitRatio[2]/100)*len(SelectedCases))
    NTrain = len(SelectedCases) - NTest - NExt
    SmiIDt, SmiID_ext = train_test_split(SmiID, test_size = NExt, random_state=42)
    SmiID_train, SmiID_test = train_test_split(SmiIDt, test_size = NTest, random_state=42)
    
    FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
    FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
    f1 = open(InputFile.split(".")[0]+"_Curated_Balanced_Training.txt", "w")
    f1.write("\t".join(FP) + "\n")
    f2 = open(InputFile.split(".")[0]+"_Curated_Balanced_Test.txt", "w")
    f2.write("\t".join(FP) + "\n")
    f3 = open(InputFile.split(".")[0]+"_Curated_Balanced_External.txt", "w")
    f3.write("\t".join(FP) + "\n")

    for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
        FP = [str(k) for k in smifp]
        FP = "\t".join(FP)
        if sid in SmiID_train:
            f1.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        if sid in SmiID_test:
            f2.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        if sid in SmiID_ext:
            f3.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
    f1.close()
    f2.close()
    f3.close()
else:
    if len(SplitRatio)==2:    
        glog.write("The Training table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_Training.txt"+ "\n")
        glog.write("The External table is presented in the file: "+  InputFile.split(".")[0]+"_Curated_Balanced_External.txt"+ "\n")
        NExt = int((SplitRatio[1]/100)*len(SelectedCases))
        SmiID_train, SmiID_ext = train_test_split(SmiID, test_size = NExt, random_state=42)
        FP = ["FP" + str(k) for k in range(len(SmilesDescription[0]))]
        FP = ["SMILEID", "SMILECOR", "CLASS"] + FP
        f1 = open(InputFile.split(".")[0]+"_Curated_Balanced_Training.txt", "w")
        f1.write("\t".join(FP) + "\n")
        f2 = open(InputFile.split(".")[0]+"_Curated_Balanced_External.txt", "w")
        f2.write("\t".join(FP) + "\n")

        for sid, smc, smicls, smifp in zip(SmileID, SmilesF, Class, SmilesDescription):
            FP = [str(k) for k in smifp]
            FP = "\t".join(FP)
            if sid in SmiID_train:
                f1.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
            if sid in SmiID_ext:
                f2.write(sid + "\t" + smc + "\t" + str(smicls) + "\t" + FP + "\n")
        f1.close()
        f2.close()
        
    else:
        print("An error ocurred during data partition, please check the Param File")
        sys.exit()

glog.write("#### END DATA PARTITION #######\n")
glog.write("############################################################################\n")
glog.write("SUCCESS COMPLETED\n")

glog.close()










