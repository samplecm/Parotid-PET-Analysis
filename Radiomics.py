import os
import numpy as np
import csv
import copy

from numpy.core.fromnumeric import sort
from Patient import Patient
import Contours
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
import Visuals


parentDirectory = os.getcwd()

def FeatureSelection(organ):
    #loop through all the patients and collect all of the radiomics and combine the data, then cluster plot

    features = np.zeros((1080,27))
    idx = 0
    for i in range(1,31):
        patient_idx = "{:02d}".format(i)
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , patient_idx)
        rightFeatureList, leftFeatureList = ReadFeaturesFromCSV(patientPath, organ)
        for j in range(18):
            features[idx,:] = rightFeatureList[j+1][1]
            features[idx+540,:] = leftFeatureList[j+1][1]
            idx += 1
    allFeatures_norm = NormalizeFeatures(CloneList(features))
    correlation_array, labels = CorrelationArray(allFeatures_norm)
    correlation_array_filtered, labels_filtered = Pairwise_Correlation_Filter(correlation_array, labels)
    #FeatureSelect_KMeans(allFeatures_norm)


    # combinedFeaturesRight_normalized = NormalizeFeatures(combinedFeaturesRight)
    # combinedFeaturesLeft_normalized = NormalizeFeatures(combinedFeaturesLeft)

    #Visuals.ClusterPlot(correlation_array, labels)
    #Visuals.ClusterPlot(correlation_array_filtered, labels_filtered)
    print("Plotted Clusters.")        
    return labels_filtered

def Get_Subseg_Features(labels, organ="parotid"):
    subseg_features = np.zeros((60, 18, len(labels))) #18 subsegs, 4 features, 30patientsx2organs
    idx = 0
    for i in range(1,31):
        patient_idx = "{:02d}".format(i)
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , patient_idx)
        rightFeatureList, leftFeatureList = ReadFeaturesFromCSV(patientPath, organ)
        j = 0
        k = 0
        while k < len(labels):
            feature_idx = 1000
            for l in range(len(rightFeatureList[0][0])):
                if labels[k] == rightFeatureList[0][0][l]:
                    feature_idx = l
                    break                    
            for j in range(len(rightFeatureList)-1):              
                subseg_features[i-1, j, k] = rightFeatureList[j+1][1][feature_idx] #5th percentile
                subseg_features[i-1+30, j, k] = leftFeatureList[j+1][1][feature_idx]
            k += 1    


    #Now condense into a statistics array
    radiomics_stats = np.zeros((2,18,len(labels))) #2 rows first for left and right parotid
    for i in range(radiomics_stats.shape[1]):
        #radiomics_stats[i,0, 0] = np.percentile(subseg_features[i,0,:], 5)  
        for k in range(len(labels)):   
            radiomics_stats[0,i,k] = np.nanmean(subseg_features[0:30, i , k])#np.percentile(subseg_features[i,0,0:30], 50)    
            radiomics_stats[1,i,k] = np.nanmean(subseg_features[30:, i , k]) 
        # #radiomics_stats[i,0, 2] = np.percentile(subseg_features[i,0,:], 95)    

        # #radiomics_stats[i,1, 0] = np.percentile(subseg_features[i,1,:], 5)       
        # radiomics_stats[0,i,1] = np.average(subseg_features[i,1,0:30]) 
        # radiomics_stats[1,i,1] = np.average(subseg_features[i,1,30:])
        # #radiomics_stats[i,1, 2] = np.percentile(subseg_features[i,1,:], 95)    

        # #radiomics_stats[i,2, 0] = np.percentile(subseg_features[i,2,:], 5)       
        # radiomics_stats[0,i,2] = np.average(subseg_features[i,2,0:30])    
        # radiomics_stats[1,i,2] = np.average(subseg_features[i,2,30:]) 
        # #radiomics_stats[i,2, 2] = np.percentile(subseg_features[i,2,:], 95)    

        # #radiomics_stats[i,3, 0] = np.percentile(subseg_features[i,3,:], 5)       
        # radiomics_stats[0,i,3] = np.average(subseg_features[i,3,0:30])   
        # radiomics_stats[1,i,3] = np.average(subseg_features[i,3,30:])  
        #radiomics_stats[i,3, 2] = np.percentile(subseg_features[i,3,:], 95)    
    return radiomics_stats    
    print("Collected population radiomics statistics for model.")
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

    array, labels = SortClusters(array)
    return array, labels              

def Pairwise_Correlation_Filter(correlation_array, labels):
    #this function removes redundant features, which are correlated over a certain threshold with another
    threshold = 0.9
    remove_list = [] #incidces of featues to remove
    num_features = correlation_array.shape[0]
    max_cor = 1
    while max_cor >= threshold:
        max_cor = 0
        max_idx = 0
        for i in range(num_features):
            if i in remove_list:
                continue
            for j in range(num_features):
                if j in remove_list or j == i:
                    continue
                cor = abs(correlation_array[i,j])
                if cor > max_cor:
                    max_cor = cor
                    max_idx = i
        if max_cor >= threshold:
            remove_list.append(max_idx)
    #now need to make a new correlation array and filter out labels that aren't used any more!
    labels_filtered = []
    for i, label in enumerate(CloneList(labels)):
        if not i in remove_list:
            labels_filtered.append(label)

    num_features -= len(remove_list)
    correlation_array_filtered = np.zeros((num_features, num_features))
    good_indices = range(correlation_array.shape[0])
    good_indices = [element for i, element in enumerate(good_indices) if i not in remove_list]

    for new_i, i in enumerate(good_indices):
        for new_j, j in enumerate(good_indices):
            correlation_array_filtered[new_i,new_j] = correlation_array[i,j]

    return correlation_array_filtered, labels_filtered        





def SortClusters(array, threshold=0.7):
    #this function sorts the correlation array to group clusters according to a certain threshold 
    #look for the feature with the most correlations, and then add all those features first, then repeat until there are no more correlations. Then search with a smaller threshold
    #until threshold is 0.2, then just add the rest.
    num_features = array.shape[0]
    sorted_features = np.zeros((num_features, num_features))

    features_moved = 0 #keep track of how many have been sorted and placed into the sorted array
    new_index_order = [] #for resorting in the end
    while features_moved < num_features:
        #first look for the row with the most correlations above the threshold
        cors = []
        most_cor = 0
        for i in range(num_features):
            if i in new_index_order:
                continue
            features = array[i,:]
            cors_indices = []
            for j in range(num_features):
                if j in new_index_order:
                    continue
                if features[j] > threshold:
                    cors_indices.append(j)

            if len(cors_indices) > most_cor:
                cors = cors_indices
                most_cor = len(cors_indices)
        #now add all these features to the first available rows in the new features array. But need to resort values first to go along with the move.
        if len(cors) < 4 and threshold > 0.6: 
            threshold -= 0.1
        elif len(cors) < 3:
            threshold -= 0.1
        new_index_order.extend(cors) 

        current_row_idx = features_moved #next available row in new array

        for i, idx in enumerate(cors):
            sorted_features[current_row_idx+ i, :]  = array[idx]    

        features_moved += most_cor

    sorted_features = Sort_Clustering_Indices(np.copy(array), new_index_order)   

    #Now get the order of feature labels after sorting. 
    filePath = "/media/caleb/WDBlue/PET_PSMA/PSMA_Analysis/SG_PETRT/03/Radiomics/rParSub_14.csv"
    with open(filePath) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            labels = row[3:]
            break
    labels = Sort_Labels(labels.copy(), new_index_order)    
    
    return sorted_features, labels 


def Sort_Clustering_Indices(array, new_index_order : list):
    #When we cluster features and reorder them, the indices swap in the 2d feature correlation array. This function performs that swapping

    sorted_array = np.zeros(array.shape)

    #go through the array and every i in cors swaps with current_row_idx + i
    for i in range(array.shape[0]):
        for j in range(array.shape[0]):
            i_idx = new_index_order[i]
            j_idx = new_index_order[j]
            sorted_array[i,j] = array[i_idx, j_idx]

    return sorted_array

    #

def Sort_Labels(labels, new_index_order):
    newLabels = []
    for i in new_index_order:
        newLabels.append(labels[i])
    return newLabels    

def NormalizeFeatures(array):
    #convert each feature to its Z value
    for i in range(27):
        avg = np.nanmean(array[:,i])
        std = np.nanstd(array[:,i])
        array[:,i] = (array[:,i] - avg)/std
    return array     

def ReadFeaturesFromCSV(patientPath, organ):
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
                    rightFeatureList[0].append(row[3:-8])
        elif left_whole_org_string in file:   
            with open(filePath) as csv_file: 
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    leftFeatureList[0].append(row[3:-8]) 

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
                        rightFeatureList[-1].append(row[3:-8])
            elif str(left_org_string+ str(s+1)+ ".csv") in file:
                with open(filePath) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        leftFeatureList[-1].append(row[3:-8])            


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



        