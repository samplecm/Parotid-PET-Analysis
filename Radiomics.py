import os
import numpy as np
import csv
import copy
from Patient import Patient
import Contours
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
import Visuals


parentDirectory = os.getcwd()

def FeatureSelection(organ):
    #loop through all the patients and collect all of the radiomics and combine the data, then cluster

    features = np.zeros((1080,35))
    idx = 0
    for i in range(1,31):
        patient_idx = "{:02d}".format(i)
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , patient_idx)
        rightFeatureList, leftFeatureList = ReadFromCSV(patientPath, organ)
        for j in range(18):
            features[idx,:] = rightFeatureList[j+1][1]
            features[idx+540,:] = leftFeatureList[j+1][1]
            idx += 1
    allFeatures_norm = NormalizeFeatures(CloneList(features))

    correlation_array = CorrelationArray(allFeatures_norm)
    #FeatureSelect_KMeans(allFeatures_norm)


    # combinedFeaturesRight_normalized = NormalizeFeatures(combinedFeaturesRight)
    # combinedFeaturesLeft_normalized = NormalizeFeatures(combinedFeaturesLeft)
    Visuals.ClusterPlot(correlation_array)
    print("Selected features")        

def CorrelationArray(features, method="spearman", load=True):
    #obtain the unsorted correlation array
    savePath = os.path.join(parentDirectory, "Statistics/correlation_array.npy")
    if load==True:
        print("Attempting to load correlation array.")  
        try:
            with open(savePath, 'rb') as f:
                array = np.load(f)
        except:
            print("Did not find saved correlation array.") 
            num_features = features.shape[1]
            array = np.zeros((num_features, num_features))
            for i in range(num_features):
                for j in range(num_features):
                    feature_list_1 = features[:,i]
                    feature_list_2 = features[:,j]
                    spearman_rank = SpearmansRankCorrelation(feature_list_1, feature_list_2)
                    array[i,j] = spearman_rank
     
    else:           
        num_features = features.shape[1]
        array = np.zeros((num_features, num_features))
        for i in range(num_features):
            for j in range(num_features):
                feature_list_1 = features[:,i]
                feature_list_2 = features[:,j]
                spearman_rank = SpearmansRankCorrelation(feature_list_1, feature_list_2)
                array[i,j] = spearman_rank

    
    with open(savePath, 'wb') as f:
        np.save(f, array)

    return array              





def NormalizeFeatures(array):
    #convert each feature to its Z value
    for i in range(35):
        avg = np.nanmean(array[:,i])
        std = np.nanstd(array[:,i])
        array[:,i] = (array[:,i] - avg)/std
    return array     

def ReadFromCSV(patientPath, organ):
    #This function grabs the radiomics features from the csv produced by DICOMautomaton's feature extractor. 
    #The data is formatted as 2 lists (one for left, one for right) containing a list for the whole ROI and a 
    # list for each subsegment, in that order. For each rois list, each feature is an element, which itself is a length-2 list, 
    #where the first element is the feature name and second element is the value. 
    csvDir = os.path.join(patientPath, "Radiomics")
    if "parotid" in organ: 
        right_org_string = "rParSub_"
        left_org_string = "lParSub_"
        left_whole_org_string = "lp_features"
        right_whole_org_string = "rp_features"
        num_subsegs = 18

    elif "subman" or "sm" in organ:
        right_org_string = "rSMSub_"
        left_org_string = "lSMSub_"
        left_whole_org_string = "ls_features"
        right_whole_org_string = "rs_features"
        num_subsegs = 8
    else:
        raise ValueError("invalid organ specified. Radiomics can be extracted for parotid or submandibular glands")
            

    leftFeatureList = [[]]
    rightFeatureList = [[]]


    #first get whole organ features
    for file in os.listdir(csvDir):
        filePath = os.path.join(csvDir, file)
        if right_whole_org_string in file:
            with open(filePath) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    rightFeatureList[0].append(row[3:])
        elif left_whole_org_string in file:   
            with open(filePath) as csv_file: 
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    leftFeatureList[0].append(row[3:]) 

#Now get subsegment features
    for s in range(num_subsegs):
        leftFeatureList.append([])
        rightFeatureList.append([])
        for file in os.listdir(csvDir):
            filePath = os.path.join(csvDir, file)
            if str(right_org_string+ str(s+1)+ ".csv") in file:
                with open(filePath) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        rightFeatureList[-1].append(row[3:])
            elif str(left_org_string+ str(s+1)+ ".csv") in file:
                with open(filePath) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        leftFeatureList[-1].append(row[3:])            


    print("Finished reading radiomics from CSVs")
    return rightFeatureList, leftFeatureList

def FeatureSelect_KMeans(features):
    kmeans = KMeans(n_clusters=5)

def TestRadiomics():
    dataset = make_blobs(n_samples=200, centers=4, n_features=2, cluster_std=1.6, random_state=50)
    points = dataset[0]
    print(points)
    #create a KMeans object
    kmeans = KMeans(n_clusters=4)
    kmeans.fit(points)

    clusters = kmeans.cluster_centers_
    print(clusters)

    y_km = kmeans.fit_predict(points)
    print(y_km)
    # plt.scatter(dataset[0][:,0], dataset[0][:,1])
    # plt.show()

    plt.scatter(points[y_km == 0,0], points[y_km == 0,1], s=50, color='red')
    plt.show()
    #print(dataset)

def CloneList(list):
    listCopy = copy.deepcopy(list)
    return listCopy    
def SpearmansRankCorrelation(list1, list2):
    #This function computes the spearmans rank coefficient for a list of suv values


    list1_Ranks = GetListRank(list1)      
    list2_Ranks = GetListRank(list2)
    n = len(list1)
    sumRankDifs = 0
    for i in range(n):
        sumRankDifs = sumRankDifs + (list1_Ranks[i] - list2_Ranks[i])**2
    rho = 1 - ((6*sumRankDifs)/(n*(n**2-1)))    
    return rho

def GetListRank(inputList):
    #This function returns a list of the same size which returns integers corresponding to the ordinal value of each item in the input list, from largest to smallest
    inputList = inputList.tolist()
    inputListOrig = CloneList(inputList)
    rankList = []
    for i in range(len(inputList)):
        rankList.append(0)
    rank = 1
    while len(inputList) > 0:
        maxVal = -1000
        maxIdx = 0
        for i in range(len(inputList)):
            if inputList[i] > maxVal:
                maxVal = inputList[i]
                maxIdx = inputListOrig.index(maxVal)
                origIdx = i
        rankList[maxIdx] = rank
        inputList.pop(origIdx)
        rank = rank + 1

    return rankList    
        
if __name__ == "__main__":
    TestRadiomics()



        