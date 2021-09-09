import csv
import os
import scipy.stats
import numpy as np 
import copy 
def SUV_Metrics():
    #this function computes the average, max and min suv inside each subsegment and whole gland of the left and parotid gland over all patients

    suvAvgs = [] #left is left, right is right parotid. thee first index is whole gland, followed by each subsegment
    suvMins = []
    suvMaxs = []
    suvSTDs = []

    for i in range(19):
        suvAvgs.append([0,0])
        suvMins.append([1e6, 1e6])
        suvMaxs.append([1e-6,1e-6])
        suvSTDs.append([0,0])

    #loop through first and get average, min and max suvs
    for i in range(1,31):
        patientPath = "/media/calebsample/Data/PET PSMA/PSMA Analysis/SG_PETRT/" + str(i)
        csvPath = os.path.join(patientPath, 'suv_stats.csv')
        with open(csvPath) as csvFile:   
            csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
            for row in range(1, 20):          
                for idx, rowData in csvData:

                    if int(idx) == int(row):
                        suvAvgs[row-1][0] = suvAvgs[row-1][0] + float(rowData[1])
                        suvAvgs[row-1][1] = suvAvgs[row-1][1] + float(rowData[2])
                        if float(rowData[1]) < suvMins[row-1][0]:
                            suvMins[row-1][0] = float(rowData[1])
                        if float(rowData[2]) < suvMins[row-1][1]:
                            suvMins[row-1][1] = float(rowData[2])    
    for row in range(1, 20):
        suvAvgs[row-1][0] = suvAvgs[row-1][0] / 30
        suvAvgs[row-1][1] = suvAvgs[row-1][1] / 30


    #now loop back through and get stds.    
    for i in range(1,31):
        patientPath = "/media/calebsample/Data/PET PSMA/PSMA Analysis/SG_PETRT/" + str(i)
        csvPath = os.path.join(patientPath, 'suv_stats.csv')
        with open(csvPath) as csvFile:   
            csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
            for row in range(1, 20):          
                for idx, rowData in csvData:

                    if int(idx) == int(row):
                        suvSTDs[row-1][0] = suvSTDs[row-1][0] + (float(rowData[1])-suvAvgs[row-1][0])**2
                        suvSTDs[row-1][1] = suvSTDs[row-1][1] + (float(rowData[2])-suvAvgs[row-1][1])**2
    for row in range(1, 20):
        suvSTDs[row-1][0] = (suvSTDs[row-1][0] / 30)**0.5
        suvSTDs[row-1][1] = (suvSTDs[row-1][1] / 30)**0.5
    


    #Now get the number of significant correlations (spearman and pearson) using a significant cutoff of p = 0.01, 0.02, 0.05, 0.1 for left and right parotids
    numSig_Spear_Left = [0,0,0,0]
    numSig_Spear_Right = [0,0,0,0]
    numSig_Pears_Left = [0,0,0,0]
    numSig_Pears_Right = [0,0,0,0]
    for i in range(1,31):
        patientPath = "/media/calebsample/Data/PET PSMA/PSMA Analysis/SG_PETRT/" + str(i)
        csvPath = os.path.join(patientPath, 'suv_stats.csv')
        
        with open(csvPath) as csvFile:   
            csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
            leftSpear, rightSpear = [float(i) for i in csvData[22][1][1:]]
            leftPears, rightPears = [float(i) for i in csvData[23][1][1:]]
            
            if leftSpear < 0.1:
                numSig_Spear_Left[0] = numSig_Spear_Left[0] + 1
            if leftSpear < 0.05: 
                numSig_Spear_Left[1] = numSig_Spear_Left[1] + 1
            if leftSpear < 0.02: 
                numSig_Spear_Left[2] = numSig_Spear_Left[2] + 1    
            if leftSpear < 0.01: 
                numSig_Spear_Left[3] = numSig_Spear_Left[3] + 1    

            if leftPears < 0.1:
                numSig_Pears_Left[0] = numSig_Pears_Left[0] + 1
            if leftPears < 0.05: 
                numSig_Pears_Left[1] = numSig_Pears_Left[1] + 1
            if leftPears < 0.02: 
                numSig_Pears_Left[2] = numSig_Pears_Left[2] + 1    
            if leftPears < 0.01: 
                numSig_Pears_Left[3] = numSig_Pears_Left[3] + 1       

            if rightSpear < 0.1:
                numSig_Spear_Right[0] = numSig_Spear_Right[0] + 1
            if rightSpear < 0.05: 
                numSig_Spear_Right[1] = numSig_Spear_Right[1] + 1
            if rightSpear < 0.02: 
                numSig_Spear_Right[2] = numSig_Spear_Right[2] + 1    
            if rightSpear < 0.01: 
                numSig_Spear_Right[3] = numSig_Spear_Right[3] + 1    

            if rightPears < 0.1:
                numSig_Pears_Right[0] = numSig_Pears_Right[0] + 1
            if rightPears < 0.05: 
                numSig_Pears_Right[1] = numSig_Pears_Right[1] + 1
            if rightPears < 0.02: 
                numSig_Pears_Right[2] = numSig_Pears_Right[2] + 1    
            if rightPears < 0.01: 
                numSig_Pears_Right[3] = numSig_Pears_Right[3] + 1        

    #now want to computer pearson and spearman coefficients using average suv values. 
    leftSubseg_avgs, rightSubseg_avgs = map(list, zip(*suvAvgs))
    leftSubseg_avgs = leftSubseg_avgs[1:]
    rightSubseg_avgs = rightSubseg_avgs[1:]
    rightSpearman = SpearmansRankCorrelation(CloneList(rightSubseg_avgs))
    leftSpearman = SpearmansRankCorrelation(CloneList(leftSubseg_avgs))
    rightSpearman_t = Get_t_value(rightSpearman, 18)
    leftSpearman_t = Get_t_value(leftSpearman, 18)
    rightSpearman_p = scipy.stats.t.sf(np.abs(rightSpearman_t), 16)
    leftSpearman_p = scipy.stats.t.sf(np.abs(leftSpearman_t), 16)
    rightPearson = PearsonsRankCorrelation(rightSubseg_avgs)
    leftPearson = PearsonsRankCorrelation(leftSubseg_avgs)
    rightPearson_t = Get_t_value(rightPearson, 18)
    leftPearson_t = Get_t_value(leftPearson, 18)
    rightPearson_p = scipy.stats.t.sf(np.abs(rightPearson_t), 16)
    leftPearson_p = scipy.stats.t.sf(np.abs(leftPearson_t), 16)
    print("Calculated SUV Stats")  
    
    #Now need to save these stats into the suv_stats csv file.
    path = "/media/calebsample/Data/PET PSMA/PSMA Analysis/PopulationStats.csv"
    with open(path, 'w') as csvFile:
        filewriter = csv.writer(csvFile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['Subsegment', 'Left-Parotid','Left-Parotid-SD', 'Right-Parotid', 'Right-Parotid-SD'])
        filewriter.writerow(["whole-gland", suvAvgs[0][0], suvSTDs[0][0], suvAvgs[0][1], suvSTDs[0][1]]) 
        for idx in range(1, len(suvAvgs)):
            filewriter.writerow([str(int(idx)), suvAvgs[idx][0], suvSTDs[idx][0], suvAvgs[idx][1], suvSTDs[idx][1]]) 
        filewriter.writerow(["Spearmans Coefficient", leftSpearman, rightSpearman])
        filewriter.writerow(["Pearsons Coefficient", leftPearson, rightPearson])
        filewriter.writerow(["Spearman Significance", leftSpearman_p, rightSpearman_p])
        filewriter.writerow(["Pearson Significance", leftPearson_p, rightPearson_p])
      
            
def Get_t_value(rho, n):
    return rho * np.sqrt((n-2)/(1-rho**2)) 
def SpearmansRankCorrelation(suvs):
    #This function computes the spearmans rank coefficient for a list of suv values

    if len(suvs) == 18: #if 18 subsegments
        importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
            1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
            0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
            0.272984146341463,   0.0484897103658537,  0.035493902439024]
    importanceRanks = GetListRank(importanceVals)      
    suvRanks = GetListRank(suvs)
    n = len(suvRanks)
    sumRankDifs = 0
    for i in range(len(suvRanks)):
        sumRankDifs = sumRankDifs + (importanceRanks[i] - suvRanks[i])**2
    rho = 1 - ((6*sumRankDifs)/(n*(n**2-1)))    
    return rho
def PearsonsRankCorrelation(suvs):
    #This function computes the spearmans rank coefficient for a list of suv values

    if len(suvs) == 18: #if 18 subsegments
        importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
            1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
            0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
            0.272984146341463,   0.0484897103658537,  0.035493902439024]
    n = len(suvs)
    sum_xy = 0
    sum_x = 0
    sum_y = 0
    sum_x2 = 0
    sum_y2 = 0

    for i in range(len(suvs)):
        sum_xy = sum_xy + suvs[i]*importanceVals[i]
        sum_x = sum_x + suvs[i]
        sum_y = sum_y + importanceVals[i]
        sum_x2 = sum_x2 + suvs[i]**2
        sum_y2 = sum_y2 + importanceVals[i]**2
    r_numerator = n*sum_xy - (sum_x * sum_y) 
    r_denominator = np.sqrt((n*sum_x2-sum_x**2)*(n*sum_y2-sum_y**2))
    r = r_numerator / r_denominator
    return r

def CloneList(list):
    listCopy = copy.deepcopy(list)
    return listCopy    
def GetListRank(inputList):
    #This function returns a list of the same size which returns integers corresponding to the ordinal value of each item in the input list, from largest to smallest
    inputListOrig = CloneList(inputList)
    rankList = []
    for i in range(len(inputList)):
        rankList.append(0)
    rank = 1
    while len(inputList) > 0:
        maxVal = 0
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

    SUV_Metrics()                
