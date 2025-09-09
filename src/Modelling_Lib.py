# -*- coding: utf-8 -*-
"""
Created on Tue Nov 2024

Simple GA-Modeling for binary outcome and smiles input.
Library of functions

@author: edu
"""

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import random
from sklearn.model_selection import KFold
import sklearn.svm
from sklearn.metrics import f1_score, confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
import pickle
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import GaussianNB

"""
#Smiles Processings and Fingerprints/Molecular descriptions
"""
def SmilesSanitaze (Smile, IgnoreSalts, IncludeIsomer, IgnoreCharge):
    if Smile.find(".")>-1:
        if IgnoreSalts:
            return None
        else:
            try:
                m = Chem.MolFromSmiles(Smile,sanitize=True)
                fragsMolAtomMappingIDs=[]
                m1 = Chem.rdmolops.GetMolFrags(m,sanitizeFrags=True, asMols=True,fragsMolAtomMapping=fragsMolAtomMappingIDs)
                lenF = [len(a) for a in fragsMolAtomMappingIDs]
                Order = np.argsort(lenF)
                MaxV = Order[-1]
                molO = m1[MaxV]
                molCharge = Chem.GetFormalCharge(molO)
                if IgnoreCharge:
                    if IncludeIsomer:
                        return Chem.MolToSmiles(molO, isomericSmiles=True, canonical=True)
                    else:
                        return Chem.MolToSmiles(molO, isomericSmiles=False, canonical=True) 
                else:
                    if molCharge==0:
                        if IncludeIsomer:
                            return Chem.MolToSmiles(molO, isomericSmiles=True, canonical=True)
                        else:
                            return Chem.MolToSmiles(molO, isomericSmiles=False, canonical=True) 
                    else:
                        return None
            except:
                return None                
    else:
        try:
            m = Chem.MolFromSmiles(Smile,sanitize=True)
            molCharge = Chem.GetFormalCharge(m)
            if IgnoreCharge:
                if IncludeIsomer:
                    return Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
                else:
                    return Chem.MolToSmiles(m, isomericSmiles=False, canonical=True)
            if molCharge==0:
                if IncludeIsomer:
                    return Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
                else:
                    return Chem.MolToSmiles(m, isomericSmiles=False, canonical=True) 
            else:
                return None
        except:
            return None
        

def ComputeDescriptors (Smile, DescriptorType, IncludeIsomer):
    try:
        m = Chem.MolFromSmiles(Smile)
        if DescriptorType==1:
            if IncludeIsomer:
                fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=True)
            else:
                fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=False)
            F = fpgen.GetFingerprint(m)
        if DescriptorType==2:
            if IncludeIsomer:
                fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
            else:
                fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)
            F = fpgen.GetFingerprint(m)
        if DescriptorType==3:
            m = Chem.MolFromSmiles(Smile)
            F = AllChem.RDKFingerprint(m)
        return F
    except:
        return None
       

"""
#Data Balancing Processings and Fingerprints/Molecular descriptions
"""

def DataBalancing (SmileDescriptors, Class, MetricCluster, Cmin):
    
    """
    #Available metrics
    # ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, 
    #‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, 
    #‘sqeuclidean’, ‘yule’, ‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’
    """
         
    X = np.array(SmileDescriptors)
    Y = np.array(Class, dtype=int)


    # BALANCING DATA

    N_1 = sum(Y==1)
    N_0 = sum(Y==0)

    if N_1 > N_0:
        Index_Include = np.where(Y==0)[0].tolist()
        x_B_index = np.where(Y==1)[0]
        x_B = X[x_B_index,:]

        print("Finding Optimal Clusters....")
        sil = []
        for k in range(Cmin, N_0):
            kmeans = KMeans(n_clusters = k).fit(x_B)
            labels = kmeans.labels_
            sil.append(silhouette_score(x_B, labels, metric = MetricCluster))
        
        K = np.array(range(Cmin, N_0))
        NC = K[np.array(sil)==max(sil)][0]
        kmeans = KMeans(n_clusters=NC, random_state=0).fit(x_B)

        TempOrder = []
        for jj in range(len(kmeans.labels_)): TempOrder.append(jj)
        TempOrder = np.array(TempOrder)
        NumberofCases = int(N_0)
        index=[]
        indexTOTAL=[]
        for i in range(NC):
            NumCases= round(NumberofCases*np.where(kmeans.labels_==i)[0].shape[0]/float(kmeans.labels_.shape[0]))
            index.append(np.random.choice(TempOrder[kmeans.labels_==i], NumCases, replace=False))
            indexTOTAL.append(TempOrder[kmeans.labels_==i])
        
        IndexT = []
        for i in index: 
            for j in i:
                IndexT.append(j)
        index = np.array(IndexT)    
        
        SelectedIndexes = index.tolist() + Index_Include
        return SelectedIndexes, NC
        
    else:
        if N_0 > N_1:
            Index_Include = np.where(Y==1)[0].tolist()
            x_B_index = np.where(Y==0)[0]
            x_B = X[x_B_index,:]
            
            sil = []
            for k in range(Cmin, N_1):
                kmeans = KMeans(n_clusters = k).fit(x_B)
                labels = kmeans.labels_
                sil.append(silhouette_score(x_B, labels, metric = MetricCluster))
            
            K = np.array(range(Cmin, N_1))
            NC = K[np.array(sil)==max(sil)][0]
            kmeans = KMeans(n_clusters=NC, random_state=0).fit(x_B)
            
            TempOrder = []
            for jj in range(len(kmeans.labels_)): TempOrder.append(jj)
            TempOrder = np.array(TempOrder)
            NumberofCases = int(N_1)
            index=[]
            indexTOTAL=[]
            for i in range(NC):
                NumCases= round(NumberofCases*np.where(kmeans.labels_==i)[0].shape[0]/float(kmeans.labels_.shape[0]))
                index.append(np.random.choice(TempOrder[kmeans.labels_==i], NumCases, replace=False))
                indexTOTAL.append(TempOrder[kmeans.labels_==i])
            
            
            IndexT = []
            for i in index: 
                for j in i:
                    IndexT.append(j)
            index = np.array(IndexT)  
            
            SelectedIndexes = index.tolist() + Index_Include
            return SelectedIndexes, NC
            
"""
#Data Reduction
"""
def SimplestDataReduction(Xtr, VarCutoff):
    VarSel = []
    for i in range(Xtr.shape[1]):
        sumT = np.sum(Xtr[:,i])
        if sumT != 0 and np.var(Xtr[:,i])>=VarCutoff:
            VarSel.append(i)
    return np.array(VarSel)    
        

"""
#Data Application Domain
"""
         

def ApplicationDomainSpaceGenerator(Xtr):
    pca = PCA()
    pca.fit_transform(Xtr)
    XPCA= pca.transform(Xtr)
    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_variance_ratio = np.cumsum(explained_variance_ratio)
    k=0
    for x,v in enumerate(cumulative_variance_ratio):
        if v>0.9:
            k=x
            break
    Ncomponents=k
    XPCA = XPCA[:,0:k]
    centroid= np.mean(XPCA, axis=0)
    centroid = centroid.reshape(1,-1)
    a=pairwise_distances(centroid, XPCA, metric='sqeuclidean')
    return pca, Ncomponents, centroid, np.max(a)

def ApplicationDomainSpaceCheck(Xte, pca, MaxDist, Ncomponents, centroid):
    XtePCA= pca.transform(Xte)
    XtePCA = XtePCA[:,0:Ncomponents]
    a=pairwise_distances(centroid, XtePCA, metric='sqeuclidean')[0]
    XDomain = []
    for i,d in enumerate(a):
        if d<=MaxDist:
            XDomain.append(i)
    return np.array(XDomain)
    
    
"""
#Data Models and Genetic Algorithms Without Ensembles
"""    
    
def pop_ini_generator(pop_size, num_features, minnv, maxnv):
    """
    Info
        pop_size        : number of individuals in the initial population
        num_features    : numeber of total variables
        minnv            : minimun number of variables
        maxnv            : maximum number of variables
        
    """
    
    Pop = np.zeros((pop_size, num_features), dtype= int)
    for i in range(pop_size):
        num_variables = random.sample(range(minnv, maxnv+1), 1)[0]
        Vector2 = np.array(random.sample(range(num_features), num_variables), dtype= int)
        Pop[i,Vector2]=1
    return Pop 



def cal_pop_fitness (pop, Ytrain, Xtrain, Ytest, Xtest, fitfunction, fitmetric):
    
    ModelFold, ACC, BCR, F1, SP, SE =[],[],[],[],[],[]
    Fitness=[]
    for curr_solution in pop:
        selected_elements_indices = np.where(curr_solution == 1)[0]
        X_Train = Xtrain[:,selected_elements_indices]
        X_Test = Xtest[:,selected_elements_indices]
    
        if fitfunction == "ML":
            model_classifier = sklearn.svm.SVC(kernel='rbf', gamma=1e-2, C=10)
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "RF":
            model_classifier = RandomForestClassifier(n_estimators=100, max_depth=4, bootstrap=False)
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "DTREE":
            model_classifier = DecisionTreeClassifier()
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "KNN":
            model_classifier = KNeighborsClassifier()
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "LDA":
            model_classifier = LinearDiscriminantAnalysis()
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "LR":
            model_classifier = LogisticRegression()
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "BNB":
            model_classifier = BernoulliNB()
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
        if fitfunction == "GNB":
            model_classifier = GaussianNB()
            model_classifier.fit(X_Train, Ytrain)
            y_pred = model_classifier.predict(X_Test)
            cm1 = confusion_matrix(Ytest,y_pred)
            
        sensitivity = cm1[0,0]/(cm1[0,0]+cm1[0,1])
        specificity = cm1[1,1]/(cm1[1,0]+cm1[1,1])
        ACC.append(accuracy_score(Ytest,y_pred))
        BCR.append((sensitivity + specificity)/2 * (1-abs(sensitivity-specificity)))
        F1.append(f1_score(Ytest,y_pred))
        SP.append(specificity)
        SE.append(sensitivity)
        ModelFold.append(model_classifier)
        
    if fitmetric == "ACC":
        Fitness = ACC
    if fitmetric == "BCR":
        Fitness = BCR
    if fitmetric == "F1":
        Fitness = F1
        
    Out = np.array([ACC, BCR, F1, SP, SE, ModelFold],dtype=object)
    Out = np.transpose(Out)
    return Fitness, Out
            

def select_mating_pool(pop, fitness, num_parents, randomselection):
    # Selecting the best individuals in the current generation as parents for producing the offspring of the next generation.   
    if num_parents % 2 == 0:
        fitnessIndex= np.argsort(fitness)
        FitnessOrdered = np.array(fitness)[fitnessIndex]
        popOrdered = pop[fitnessIndex,:]
        weight_norm = FitnessOrdered / FitnessOrdered.sum()
        weight_comu = weight_norm.cumsum()
        selectedParents = []
        if randomselection:
            k=0
            while k<num_parents:
                a = random.choices(range(pop.shape[0]), k=1)[0]
                if a not in selectedParents:
                    selectedParents.append(a)
                    k=k+1            
        else:
            k=0
            while k<num_parents:
                a = random.choices(range(pop.shape[0]), cum_weights=weight_comu, k=1)[0]
                if a not in selectedParents:
                    selectedParents.append(a)
                    k=k+1

        selectedParents = np.array(selectedParents)
        parents= popOrdered[selectedParents, :]

        return parents
    else:
        return "ERROR"
   
        
def crossover(parents, offspring_size, max_features):
    #Single point crossover 2 child from 2 parents
    offspring=[]
    if parents.shape[0]==2:
        Pass = False
        while Pass == False:
            point = random.sample(range(1,parents.shape[1]), 1)[0]
            P_1 , P_2 = parents[0,:].tolist() , parents[1,:].tolist()
            child_1 = np.concatenate((P_1[:point],P_2[point:]))
            child_2 = np.concatenate((P_2[:point],P_1[point:]))
            if np.sum(child_1) <= max_features and  np.sum(child_2) <= max_features:
                offspring.append(child_1)
                offspring.append(child_2)
                Pass = True
    else:        
        for i in range(offspring_size//2):
            Pass = False
            while Pass == False:
                PairParents = np.random.choice(range(parents.shape[0]), 2, replace=False)
                point = random.sample(range(1,parents.shape[1]), 1)[0]
                P_1 , P_2 = parents[PairParents[0],:].tolist() , parents[PairParents[1],:].tolist()
                child_1 = np.concatenate((P_1[:point],P_2[point:]))
                child_2 = np.concatenate((P_2[:point],P_1[point:]))
                if np.sum(child_1) <= max_features and  np.sum(child_2) <= max_features:
                    offspring.append(child_1)
                    offspring.append(child_2)
                    Pass = True    
    return offspring
        

def mutation(offspring_crossover, prob_mutation, max_features):
    mutated_offspring = offspring_crossover.copy()
    mutation_number = int(len(offspring_crossover)*len(offspring_crossover[0])*prob_mutation)
    for k in range(mutation_number):
        i = np.random.randint(0,len(mutated_offspring))
        j = np.random.randint(0,len(mutated_offspring[0]))
        if mutated_offspring[i][j]==0 and np.sum(mutated_offspring[i]) < max_features:
            mutated_offspring[i][j]=1
        if mutated_offspring[i][j]==1 and np.sum(mutated_offspring[i]) > 2:
            mutated_offspring[i][j]=0
    return np.array(mutated_offspring)

"""
Ensembles functions
"""

def ini_models_generator(N_Initial_Models, InitialACC, min_features_ensemble, 
                         max_features_ensemble, X_train, X_test, y_train, y_test, FitnessMetric):
    
    num_features = X_train.shape[1]
    ModelFold, ACC, BCR, F1, SP, SE =[],[],[],[],[],[]

    Pop=[]
    i=0
    while i < N_Initial_Models:
        num_variables = random.sample(range(min_features_ensemble, max_features_ensemble+1), 1)[0]
        Vector2 = np.array(random.sample(range(num_features), num_variables), dtype= int)
        Xtrain = X_train[:, Vector2]
        Xtest = X_test[:, Vector2]
        
        k=np.random.randint(0,4)
        
        # if k ==0:
        #     model_classifier = sklearn.svm.SVC(kernel='rbf', probability=True, gamma=1e-2, C=10)
        # elif k==1:
        #     model_classifier = RandomForestClassifier(n_estimators=100, max_depth=4, bootstrap=False)
        # elif k==2:
        #     model_classifier = DecisionTreeClassifier()
        # else:
        #     model_classifier = KNeighborsClassifier()
            
        if k==0:
            model_classifier = RandomForestClassifier(n_estimators=100, max_depth=4, bootstrap=False)
        elif k==1:
            model_classifier = DecisionTreeClassifier()
        elif k==2:
            model_classifier = BernoulliNB()
        elif k==3:
            model_classifier = GaussianNB()
        else:
            model_classifier = KNeighborsClassifier()
        
        model_classifier.fit(Xtrain, y_train)
        y_pred = model_classifier.predict(Xtest)
        cm1 = confusion_matrix(y_test,y_pred)
        sensitivity = cm1[0,0]/(cm1[0,0]+cm1[0,1])
        specificity = cm1[1,1]/(cm1[1,0]+cm1[1,1])
        acu = accuracy_score(y_test,y_pred)
        bcrt = (sensitivity + specificity)/2 * (1-abs(sensitivity-specificity))
        
        if FitnessMetric == "ACC":
            aa=acu
        if FitnessMetric == "BCR":
            aa=bcrt
        if FitnessMetric == "F1":
            aa=f1_score(y_test,y_pred)

        if aa>=InitialACC:
            ACC.append(accuracy_score(y_test,y_pred))
            BCR.append((sensitivity + specificity)/2 * (1-abs(sensitivity-specificity)))
            F1.append(f1_score(y_test,y_pred))
            SP.append(specificity)
            SE.append(sensitivity)
            ModelFold.append(model_classifier)
            Pop.append(Vector2)
            print(i)
            i=i+1

    return Pop, ACC, BCR, F1, SP, SE, ModelFold


def cal_pop_fitness_ensemble (pop, X_test, y_test, ModelFold, ModelsVariables, FitnessMetric):
    
    ACC, BCR, F1, SP, SE =[],[],[],[],[]
    Fitness=[]
    for curr_solution in pop:
        selected_elements_indices = np.where(curr_solution == 1)[0]
        ModelsSelected = [ModelFold[i] for i in selected_elements_indices]
        ModelsVariablesSelected = [ModelsVariables[i] for i in selected_elements_indices] 
        YpAll = []
        for model, variable in zip(ModelsSelected, ModelsVariablesSelected):
            Xt = X_test[:, variable]
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
        cm1 = confusion_matrix(y_test,YpBin)
        sensitivity = cm1[0,0]/(cm1[0,0]+cm1[0,1])
        specificity = cm1[1,1]/(cm1[1,0]+cm1[1,1])
        ACC.append(accuracy_score(y_test,YpBin))
        BCR.append((sensitivity + specificity)/2 * (1-abs(sensitivity-specificity)))
        F1.append(f1_score(y_test,YpBin))
        SP.append(specificity)
        SE.append(sensitivity)
        
        if FitnessMetric == "ACC":
            aa=accuracy_score(y_test,YpBin)
        if FitnessMetric == "BCR":
            aa=(sensitivity + specificity)/2 * (1-abs(sensitivity-specificity))
        if FitnessMetric == "F1":
            aa=f1_score(y_test,YpBin)
            
        Fitness.append(aa)
    Out = np.array([ACC, BCR, F1, SP, SE],dtype=object)
    Out = np.transpose(Out)
        
    return Fitness, Out  