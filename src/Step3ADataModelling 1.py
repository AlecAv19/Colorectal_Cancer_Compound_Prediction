# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 07:41:08 2025

@author: edu

This code is part of the Pipeline-Codes
Part 3A. Data Modelling using the Split of three Files, Training, Test and External

Input: 
    Type : Only the param file

Read carefully the ParamFile (Step2_Parameters.txt)    
    
Processing Format:
    python Step1SmilesBalancingPartition.py -p ParamFile
    python Step1SmilesBalancingPartition.py -h
"""

from Modelling_Lib import  ApplicationDomainSpaceGenerator, ApplicationDomainSpaceCheck, pop_ini_generator, cal_pop_fitness, select_mating_pool, crossover, mutation, ini_models_generator,cal_pop_fitness_ensemble
import argparse
import sys
import numpy as np
import pickle
from sklearn.metrics import f1_score, confusion_matrix
from sklearn.metrics import accuracy_score


def ModifyPopulation (new_population, fitness, Others, offspring_mutation, fitness2, Others2):
    
    PopulationF1 = fitness.copy()
    Others1 = Others.copy()
    Population1 = new_population.copy()
    
    OrderF2 = np.argsort(np.array(fitness2))[::-1]
    PopulationF2 = np.array(fitness2)[OrderF2].tolist()
    Others2 = Others2[OrderF2, :]
    Population2 = offspring_mutation[OrderF2, :]
    
    for i,_ in enumerate(PopulationF2):
        FDiff = []
        for j,_ in enumerate (PopulationF1):
            FDiff.append(PopulationF2[i]-PopulationF1[j])
        FDiffmax = max(FDiff)
        if FDiffmax>0:
            k = FDiff.index(FDiffmax)
            PopulationF1[k] = PopulationF2[i]
            Others1[k,:] = Others2[i,:]
            Population1[k,:] = Population2[i,:]
            
    return Population1, PopulationF1, Others1

def ReadingParamFile (Filename):
    f = open(Filename, "r")
    FinalListElements = 0
    for line in f:
        if line != "\n" and line[0] != "#":
            Temp = line.split("\n")[0].split(" = ")
            if Temp[0] == "Training_File":
                if Temp[1]!="":
                    Training_File=Temp[1]
                    FinalListElements=FinalListElements+1
                else:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "Test_File":
                if Temp[1]!="":
                    Test_File=Temp[1]
                    FinalListElements=FinalListElements+1
                else:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "External_File":
                if Temp[1]!="":
                    External_File=Temp[1]
                    FinalListElements=FinalListElements+1
                else:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "Normalization_File":
                if Temp[1]!="":
                    Normalization_File=Temp[1]
                    FinalListElements=FinalListElements+1
                else:
                    Normalization_File=None
            if Temp[0] == "Fitness_Functions":
                if Temp[1]!="":
                    Fitness_Functions=Temp[1].split(", ")
                    FinalListElements=FinalListElements+1
                else:
                    Fitness_Functions=None
            if Temp[0] == "Fitness_Metric":
                if Temp[1]!="":
                    Fitness_Metric=Temp[1]
                    FinalListElements=FinalListElements+1
                else:
                    Fitness_Metric=None
            if Temp[0] == "population_size":
                try:
                    population_size=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "min_features":
                try:
                    min_features=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "max_features":
                try:
                    max_features=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "num_generations":
                try:
                    num_generations=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "num_parents_mating":
                try:
                    num_parents_mating=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "num_children":
                try:
                    num_children=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "prob_mutation":
                try:
                    prob_mutation=float(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "N_Initial_Models":
                try:
                    N_Initial_Models=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "Initial_metrics":
                try:
                    Initial_metrics=float(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "min_features_ensemble":
                try:
                    min_features_ensemble=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "max_features_ensemble":
                try:
                    max_features_ensemble=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "population_ensemble_size":
                try:
                    population_ensemble_size=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "min_models_ensemble":
                try:
                    min_models_ensemble=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "max_models_ensemble":
                try:
                    max_models_ensemble=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None
            if Temp[0] == "num_generation_ensembles":
                try:
                    num_generation_ensembles=int(Temp[1])
                    FinalListElements=FinalListElements+1
                except:
                    print("Check the correct spelling (i.e. spaces and ending) in " + Temp[0])
                    return None

    f.close()
    return Training_File, Test_File, External_File, Normalization_File, Fitness_Functions, Fitness_Metric, population_size, min_features, max_features, num_generations, num_parents_mating, num_children, prob_mutation, N_Initial_Models,Initial_metrics, min_features_ensemble,  max_features_ensemble, population_ensemble_size, min_models_ensemble, max_models_ensemble, num_generation_ensembles 


parser = argparse.ArgumentParser(description="Step 3A of the Pipeline: Data Modelling using 3 Splits of the Data")
parser.add_argument("-p", type=str, help="Input file with computation parameters.")

if len(sys.argv) == 1:  # No arguments provided
    sys.argv = ['Step3A.py', '-p', 'Step3AParamFile.txt']
    
args = parser.parse_args()

if args.p == None:
    print("Missing Parameter file")
    sys.exit()

#Reading Parameter File
ParamFile = args.p
print(ReadingParamFile (ParamFile))
Training_File, Test_File, External_File, Normalization_File, Fitness_Functions, Fitness_Metric, population_size, min_features, max_features, num_generations, num_parents_mating, num_children, prob_mutation, N_Initial_Models,Initial_metrics, min_features_ensemble,  max_features_ensemble, population_ensemble_size, min_models_ensemble, max_models_ensemble, num_generation_ensembles = ReadingParamFile (ParamFile)

#Reading Training File
print("Loading Training File")

f = open(Training_File, "r")
glog = open(Training_File.split(".")[0]+"_Step3A_logFile.log", "w")
X, Y=[], []
N = 0
for line in f:
    if N > 0:
        Temp = line.split("\n")[0].split("\t")
        Xtemp = Temp[3:]
        Y.append(int(Temp[2]))
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
X_train = np.array(X)
Y_train = np.array(Y)

#Reading Test File
print("Loading Test File")
f = open(Test_File, "r")
TestIDs = []
X, Y=[], []
N = 0
for line in f:
    if N > 0:
        Temp = line.split("\n")[0].split("\t")
        TestIDs.append(Temp[0])
        Y.append(int(Temp[2]))
        Xtemp = Temp[3:]
        TypeTemp = " ".join(Xtemp)
        if TypeTemp.find(".")>-1:
            Xtemp = [float(a) for a in Xtemp]
        else:
            Xtemp = [int(a) for a in Xtemp]
        X.append(Xtemp)

    N=N+1
f.close()
X_test = np.array(X)
Y_test = np.array(Y)

#Reading External File
print("Loading External File")
f = open(External_File, "r")
ExtIDs = []
X, Y=[], []
N = 0
for line in f:
    if N > 0:
        Temp = line.split("\n")[0].split("\t")
        ExtIDs.append(Temp[0])
        Xtemp = Temp[3:]
        Y.append(int(Temp[2]))
        TypeTemp = " ".join(Xtemp)
        if TypeTemp.find(".")>-1:
            Xtemp = [float(a) for a in Xtemp]
        else:
            Xtemp = [int(a) for a in Xtemp]
        X.append(Xtemp)

    N=N+1
f.close()
X_ext = np.array(X)
Y_ext = np.array(Y)

#Reading Normalization File
if Normalization_File != None:
    print("Loading Normalization File and processing")
    with open(Normalization_File, 'rb') as f:
        NormalizationData = pickle.load(f)
    Scaler = NormalizationData[0]
    FP_Selection = NormalizationData[1]
    if Scaler != None:
        X_train = Scaler.transform(X_train)
        X_test = Scaler.transform(X_test)
        X_ext = Scaler.transform(X_ext)
    if FP_Selection != None:
        Selected = [i for i, fp in enumerate(FPAll) if fp in FP_Selection]
        Selected = np.array(Selected)
        X_train = X_train[:, Selected]
        X_test = X_test[:, Selected]
        X_ext = X_ext[:, Selected]

#Application Domain
glog.write("############################################################################\n")
glog.write("#### BEGIN STEP 3A DATA MODELLING #######\n")
glog.write("#### Aplication Domain using PCA of the Normalized and Selected Variables #######\n")
print("Applying Application Domain")
pca, Ncomponents, centroid, MaxDist = ApplicationDomainSpaceGenerator(X_train)
glog.write("#### Number of Identified Principal Componets: "+ str(Ncomponents) +" #######\n")
glog.write("#### Maximal Distance to Centroid: "+ str(MaxDist) +" #######\n")
glog.write("#### Printing the Centroid Coordinates #######\n")
for a in centroid:
    glog.write(str(a) +"\n")
glog.write("\n")
glog.write("#### All parameters of the Application Domain are saved in : " + "Step3A_ApplicationDomain_Conditions.pkl" + " #######\n")
ADSAved = [pca, MaxDist, Ncomponents, centroid]
with open("Step3A_ApplicationDomain_Conditions.pkl",'wb') as f:
    pickle.dump(ADSAved,f)
        

X_Test_Domain = ApplicationDomainSpaceCheck(X_test, pca, MaxDist, Ncomponents, centroid)
X_Ext_Domain = ApplicationDomainSpaceCheck(X_ext, pca, MaxDist, Ncomponents, centroid)

glog.write("#### After applying the Application Domain to the Test and External the conclusions are: #######\n")
glog.write("Aplication Domain in Test: " + str(len(X_Test_Domain)) + " of " + str(X_test.shape[0]) +"\n")
glog.write("Aplication Domain in External: " + str(len(X_Ext_Domain)) + " of " + str(X_ext.shape[0]) +"\n")
print("Aplication Domain in Test: " + str(len(X_Test_Domain)) + " of " + str(X_test.shape[0]))
print("Aplication Domain in External: " + str(len(X_Ext_Domain)) + " of " + str(X_ext.shape[0]))

glog.write("#### Printing the IDs of the FInal Selected in Test after Application Domain #######\n")
for a in X_Test_Domain:
    glog.write(TestIDs[a] +"\n")
glog.write("\n")
glog.write("#### Printing the IDs of the FInal Selected in External after Application Domain #######\n")
for a in X_Ext_Domain:
    glog.write(ExtIDs[a] +"\n")
glog.write("\n")

X_test = X_test[X_Test_Domain,:]
X_ext = X_ext[X_Ext_Domain,:]
Y_test = Y_test[X_Test_Domain]
Y_ext = Y_ext[X_Ext_Domain]

#Genetic Algorithms and Serching
print("Starting Genetic Algorithms....")

for FitnessFunction in Fitness_Functions:
    print("Exploring " + FitnessFunction)
    num_feature = X_train.shape[1]
    new_population = pop_ini_generator(population_size, num_feature, min_features, max_features)
    fitness, Others = cal_pop_fitness(new_population, Y_train, X_train, Y_test, X_test, FitnessFunction, Fitness_Metric)
    
    i=1
    for generation in range(num_generations):
        k= fitness.index(np.max(fitness))
        print("Generation : ", i, "Best results : ", fitness[k], np.sum(new_population[k]), Others[k, 0])
        ContinueV = False
        offspring_mutation=[]
        while not ContinueV:
            parents = select_mating_pool(new_population, fitness, num_parents_mating, True)
            offspring_crossover = crossover(parents, num_children, max_features)
            offspring_mutationT = mutation(offspring_crossover, prob_mutation, max_features)
            for offP in offspring_mutationT:
                bb= [a for a in new_population if np.sum(a==offP)==num_feature]
                if bb ==[]:
                    offspring_mutation.append(offP)
                    ContinueV=True
        offspring_mutation=np.array(offspring_mutation)
        fitness2, Others2 = cal_pop_fitness(offspring_mutation, Y_train, X_train, Y_test, X_test, FitnessFunction, Fitness_Metric)
        new_population,fitness,Others = ModifyPopulation (new_population, fitness, Others, offspring_mutation, fitness2, Others2)
        i=i+1

    print("Evaluating the External Data")
    PredExternal = []
    for a, b in zip(Others, new_population):
        selected_elements_indices = np.where(b == 1)[0]
        yp = a[5].predict(X_ext[:,selected_elements_indices])
        cm1 = confusion_matrix(Y_ext,yp)
        sensitivity = cm1[0,0]/(cm1[0,0]+cm1[0,1])
        specificity = cm1[1,1]/(cm1[1,0]+cm1[1,1])
        PredExternal.append([accuracy_score(Y_ext,yp), (sensitivity + specificity)/2 * (1-abs(sensitivity-specificity)), f1_score(Y_ext,yp), specificity,sensitivity ])

    OutputFileGlobal = "Step3A_GA_Results_" + FitnessFunction + "_" + Fitness_Metric + ".txt"
    f = open(OutputFileGlobal, "w")
    f.write("GA_ACC" + "\t" + "GA_BCR" + "\t""GA_F1" + "\t""GA_SP" + "\t""GA_SE" + "\t" + 
            "EXT_ACC" + "\t" + "EXT_BCR" + "\t""EXT_F1" + "\t" + "EXT_SP" + "\t" + "EXT_SE" + "\n")
    PredMierda = [[str(b[0]), str(b[1]),str(b[2]),str(b[3]),str(b[4]), str(a[0]), str(a[1]), str(a[2]), str(a[3]), str(a[4])] for b, a in zip(Others, PredExternal)]   
    for line in PredMierda:
        f.write("\t".join(line)+"\n")
    f.close()
    
    glog.write("#### Data is saved in file: "+OutputFileGlobal+" #######\n")

    #Saving Models and Data
    InitialParameters = [OutputFileGlobal, population_size,
                         min_features, max_features, num_generations, num_parents_mating,
                         num_children, prob_mutation, FitnessFunction, Fitness_Metric]
    if Normalization_File != None:
        DataModeling = [Scaler, Selected, FP_Selection, pca, Ncomponents, centroid, MaxDist, X_train, X_test, X_ext, Y_train, Y_test, Y_ext]
    else:
        DataModeling = [pca, Ncomponents, centroid, MaxDist, X_train, X_test, X_ext, Y_train, Y_test, Y_ext]

    DataGA = [np.array(fitness), Others, new_population, np.array(PredExternal)]
    
    DataSave = [InitialParameters, DataModeling, DataGA]
    OutputFileGlobal = "Step3A_GA_Results_" + FitnessFunction + "_" + Fitness_Metric + ".pkl"
    with open(OutputFileGlobal,'wb') as f:
        pickle.dump(DataSave,f)

    glog.write("#### Data is saved in file: "+OutputFileGlobal+" #######\n")


#Ensemble Modelling

glog.write("#### Proceding to Ensembles models #######\n")
print("Starting the Ensemble Modelling")


ModelsVariables, ACC, BCR, F1, SP, SE, ModelsBest = ini_models_generator(N_Initial_Models, Initial_metrics, min_features_ensemble, 
                         max_features_ensemble, X_train, X_test, Y_train, Y_test, Fitness_Metric)

new_population = pop_ini_generator(population_ensemble_size, N_Initial_Models, min_models_ensemble, max_models_ensemble)
fitness, Others = cal_pop_fitness_ensemble (new_population, X_test, Y_test, ModelsBest, ModelsVariables, Fitness_Metric)
i=1
for generation in range(num_generation_ensembles):
    k= fitness.index(np.max(fitness))
    print("Generation : ", i, "Best results : ", fitness[k], np.mean(fitness), Others[k, 0])
    #print(TestTest(fitness, new_population, ModelsBest, ModelsVariables, X_ext, y_ext))
    ContinueV, ContinueF = False, False
    offspring_mutation=[]
    while not ContinueV and not ContinueF:        
        parents = select_mating_pool(new_population, fitness, num_parents_mating, False)
        offspring_crossover = crossover(parents, num_children, max_models_ensemble)
        offspring_mutationT = mutation(offspring_crossover, prob_mutation, max_models_ensemble)
        for offP in offspring_mutationT:
            bb= [a for a in new_population if np.sum(a==offP)==N_Initial_Models]
            if bb ==[]:
                offspring_mutation.append(offP)
                ContinueV=True
        if offspring_mutation != []:
            offspring_mutation=np.array(offspring_mutation)
            fitness2, Others2 = cal_pop_fitness_ensemble (offspring_mutation, X_test, Y_test, ModelsBest, ModelsVariables, Fitness_Metric)
            if np.max(fitness2) > np.min(fitness):
                ContinueF=True
    new_population,fitness,Others = ModifyPopulation (new_population, fitness, Others, offspring_mutation, fitness2, Others2)
    i=i+1

#Models Evaluation External
PredExternal = []
for b in new_population:
    selected_elements_indices = np.where(b == 1)[0]
    
    ModelsSelected = [ModelsBest[i] for i in selected_elements_indices]
    ModelsVariablesSelected = [ModelsVariables[i] for i in selected_elements_indices]
    
    YpAll = []
    for model, variable in zip(ModelsSelected, ModelsVariablesSelected):
        Xt = X_ext[:, variable]
        yP = model.predict_proba(Xt)[:,1]
        YpAll.append(yP)
    YpAll=np.array(YpAll)
    YpAll=np.transpose(YpAll)
    YpAllm=np.mean(YpAll, axis=1)
    YpBin=[0]*len(YpAllm)
    for i,a in enumerate(YpAllm):
        if a>0.5:
            YpBin[i]=1
    YpBin = np.array(YpBin, dtype=int)
    cm1 = confusion_matrix(Y_ext,YpBin)
    sensitivity = cm1[0,0]/(cm1[0,0]+cm1[0,1])
    specificity = cm1[1,1]/(cm1[1,0]+cm1[1,1])
    PredExternal.append([accuracy_score(Y_ext,YpBin), (sensitivity + specificity)/2 * (1-abs(sensitivity-specificity)), f1_score(Y_ext,YpBin), specificity,sensitivity ])

OutputFileGlobalEnsemble = "Step3A_Ensembles_GA" + "_" + Fitness_Metric + ".txt"
f = open(OutputFileGlobalEnsemble, "w")
f.write("GA_ACC" + "\t" + "GA_BCR" + "\t""GA_F1" + "\t""GA_SP" + "\t""GA_SE" + "\t" + 
        "EXT_ACC" + "\t" + "EXT_BCR" + "\t""EXT_F1" + "\t""EXT_SP" + "\t" + "EXT_SE" + "\n")
PredMierda = [[str(b[0]), str(b[1]),str(b[2]),str(b[3]),str(b[4]), str(a[0]), str(a[1]), str(a[2]), str(a[3]), str(a[4])] for b, a in zip(Others, PredExternal)]   
for line in PredMierda:
    f.write("\t".join(line)+"\n")
f.close()
glog.write("#### Data is saved in file: "+OutputFileGlobal+" #######\n")


#Saving Models and Data

InitialParameters = [OutputFileGlobalEnsemble,N_Initial_Models,
                     population_ensemble_size,Initial_metrics,min_features_ensemble,max_features_ensemble,
                     min_models_ensemble, max_models_ensemble, num_generation_ensembles]

if Normalization_File != None:
    DataModeling = [Scaler, Selected, FP_Selection, pca, Ncomponents, centroid, MaxDist, X_train, X_test, X_ext, Y_train, Y_test, Y_ext]
else:
    DataModeling = [pca, Ncomponents, centroid, MaxDist, X_train, X_test, X_ext, Y_train, Y_test, Y_ext]

DataGA = [np.array(fitness), Others, new_population, np.array(PredExternal)]

OutputFileGlobalEnsemble = "Step3A_Ensembles_GA" + "_" + Fitness_Metric + ".pkl"
DataSave = [InitialParameters, DataModeling, DataGA, ModelsBest, ModelsVariables]
with open(OutputFileGlobalEnsemble,'wb') as f:
    pickle.dump(DataSave,f)





        
        
glog.write("#### END STEP 3A DATA MODELLING #######\n")
glog.write("############################################################################\n")
glog.write("SUCCESS COMPLETED\n")

glog.close()








