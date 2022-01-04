import csv
from datetime import date
import os
from scipy.sparse import data
import scipy.stats
import numpy as np 
import copy
from NormalizeImportance import SubImportance 
from Patient import Patient
from Contours import Contours
import GetImageData
import pickle 
from rtstruct_builder import RTStructBuilder
from rtstruct import RTStruct
from rtutils import ROIData
import pydicom
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from numpy import mean
from numpy import absolute
from numpy import sqrt
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import itertools
parentDirectory = os.getcwd()

def Model_Selection(radiomics_data, labels):
    #Now need to test all combinations of 2,3,4,5, features for creating a model with a degree of 1,2,3,4,5
    #keep track of kfold loss and sort, and save in spreadsheet.
    degrees = np.linspace(1, 5, 5, dtype=int)
    num_features = len(labels)
    features = np.linspace(0, num_features-1, num_features)
    radiomics_data = np.round(radiomics_data, 4)
    y = []
    params_list = []
    for num_predictors in range(2, 6):
        feature_combos = itertools.combinations(features, num_predictors)
        for degree in degrees:
            for feature_indices in feature_combos:
                score = KFold_Validation(radiomics_data, feature_indices, k=9, degree=degree)
                params = [degree, [labels[int(feature_idx)] for feature_idx in feature_indices], score]
                params_list.append(params)
                y.append(score)
    params_list = sorted(params_list, key=lambda x: x[2])
    with open(os.path.join(parentDirectory, "model_selection.csv"), 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["Degree", "Features", "9-fold Loss"])
        csv_writer.writerows(params_list)


class Model:

    def __init__(self, radiomics_data, degree=2):
        self.radiomics_data = radiomics_data
        self.degree = degree
        self.fit(degree)
        
    def fit(self, subsegs="all"):
        radiomics_stats = self.radiomics_data
            #This function creates a predictive model for subsegment importance using radiomics features.
        importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
                    1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
                    0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
                    0.272984146341463,   0.0484897103658537,  0.035493902439024]
        points = []
        for i in range(radiomics_stats.shape[1]):
            points.append([radiomics_stats[0, i, 0], radiomics_stats[0, i, 1], radiomics_stats[0, i, 2], radiomics_stats[0, i, 3]])    
        for i in range(radiomics_stats.shape[1]):
            points.append([radiomics_stats[1, i, 0], radiomics_stats[1, i, 1], radiomics_stats[1, i, 2], radiomics_stats[1, i, 3]])   
        poly = PolynomialFeatures(self.degree)     
        #transform the data for fitting  
        points_transformed = poly.fit_transform(points)

        model = linear_model.LinearRegression()
        y=importanceVals*2
        model.fit(points_transformed, y)
        self.predictor = model

    def predict(self, radiomics_features):
        poly = PolynomialFeatures(self.degree)
        points = np.array(radiomics_features)
        if len(points.shape) == 1:
            points = points.reshape(1,-1)
        return self.predictor.predict(poly.fit_transform(points))


    #testSubseg = model.predict(poly.fit_transform(np.array(points[15]).reshape(1,-1)))


def KFold_Validation(radiomics_data, feature_indices, k=9, degree=2):
    #takes population radiomics data and performs a k-fold cross validation   
    importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
                    1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
                    0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
                    0.272984146341463,   0.0484897103658537,  0.035493902439024]
    y = importanceVals * 2             
    radiomics_data = np.reshape(radiomics_data, (36, 11))
    data_dict = {}
    data_cols = []
    data_dict['y'] = y
    for i, feature_i in enumerate(feature_indices):
        if np.isnan(np.min(radiomics_data[:, int(feature_i)])):
            print("Error, found nan in array for feature " + str(feature_i))
            return 1000
        data_dict[str('x' + str(i))] = radiomics_data[:, int(feature_i)]
        data_cols.append(str('x' + str(i)))
    
    data_frame = pd.DataFrame(data_dict)
    # x = data_frame[['x1', 'x2', 'x3', 'x4']]

    x = data_frame[data_cols]
    y = data_frame['y']
    poly = PolynomialFeatures(degree)     
    #transform the data for fitting  
    points_transformed = poly.fit_transform(x)

    cv = KFold(n_splits=9, random_state=1, shuffle=True)
    model = LinearRegression()
    model.fit(points_transformed, y)
    scores = cross_val_score(model, points_transformed, y, scoring="neg_mean_absolute_error", cv=cv, n_jobs=1)

    # x_train, x_test, y_train, y_test = train_test_split(np.reshape(np.array([data_frame['x1'], data_frame['x2'], data_frame['x3'], data_frame['x4']]), (36, 4)),y,train_size=(8/9))
    #print("KFold MAE: " + str(np.mean(np.abs(scores))))
    return np.mean(np.abs(scores))

    # #Loop through degrees 1 to 3 to compare performance
    # for d in range(1, 3):






def DicomSaver(patientPath, organs : list):

    contoursList = []
    with open(os.path.join(patientPath, "PatientData.txt"), "rb") as fp:
        patient : Patient = pickle.load(fp)
    CTArray = patient.CTArray
    num_CTSlices = CTArray.shape[1]    
    for string in organs:    
        if "parotid" in string.lower():
            leftPar : Contours = patient.LeftParotid
            leftParSubs = leftPar.segmentedContours18
            leftParName = "LPar"

            rightPar : Contours = patient.RightParotid
            rightParSubs = rightPar.segmentedContours18
            rightParName = "RPar"

            contoursList.append([leftParSubs, leftParName])
            contoursList.append([rightParSubs, rightParName])
        if "submandibular" in string.lower():
            leftSub : Contours = patient.LeftSubmandibular
            leftSubSubs = leftSub.segmentedContours8
            leftSubName = "LSub"

            rightSub : Contours = patient.RightSubmandibular
            rightSubSubs = rightSub.segmentedContours8
            rightSubName = "RSub"

            contoursList.append([leftSubSubs, leftSubName])
            contoursList.append([rightSubSubs, rightSubName])    

    structFile = None 

    for fileName in patientPath:
        if "Subsegs" in fileName:
            structFile = fileName  
    patientFolders = os.listdir(patientPath)
    for folder in patientFolders:
        if "CT_" in folder and "STRUCT" not in folder:
            seriesPath = os.path.join(patientPath, folder)
            #referenceCT_Meta = pydicom.dcmread(os.listdir(os.path.join(patientPath, folder)[0])

    #create a new struct file if there wasn't one provided
    if structFile is None:
        #first get the reference series from the pet folder: 
        
        rtStruct = RTStructBuilder.create_new(dicom_series_path = seriesPath)
        newStructFile = "Subsegs"
    else:
        structPath = os.path.join(patientPath, structFile)
        #load existing RT Struct
        rtStruct = RTStructBuilder.create_from(dicom_series_path = patientPath, rt_struct_path = structPath)
        newStructFile = structFile.split(".dcm")[0]

    #create a list of colors for the contours 
    colorList= [
        [255, 0, 255], #magenta
        [0, 235, 235], #teal
        [255, 255, 0], #yellow
        [255, 0, 0], #red
        [255, 175, 0], #orange
        [160, 32, 240], #purple
        [255, 140, 190], #pink
        [0, 0, 255], #blue
    ]

    for i, item in enumerate(contoursList):
        
        contours = item[0]
        #reshape contours so that it is compatible with dicom files
        for j, subsegment in enumerate(contours):
            roiName = "subseg_" + item[1].replace(" ", "")
            roiName = roiName[0:14] + str(j+1)
            newContours = []    
            for i in range(num_CTSlices):
                z_slice = CTArray[3,i,0,0]
                newPoints = []
                for slice in subsegment:
                    if len(slice) == 0 :
                        continue
                    if abs(int(round(slice[0][2], 2)*100) - int(round(z_slice, 2)*100)) < 2:
                        for point in slice:
                            newPoints.append(round(point[0],2))
                            newPoints.append(round(point[1],2))
                            newPoints.append(round(point[2],2))
                newContours.append(newPoints)    

            #assign the contour a color

            organColor = colorList[j % len(colorList)]


            #add the ROI
            rtStruct.add_roi(contours = newContours, color = organColor, name = roiName)                                

    #save the ROI to a new struct file
    rtStruct.save(str(os.path.join(patientPath, newStructFile)))
    print("Saved structure file.")

def SubmandibularSUVAnalysis(patient : Patient):   
    #this function takes a pet suv array and a structure and computes the mean suv within the structure    
    
    pet_array = patient.PETArray 
    lSub = patient.LeftSubmandibular
    rSub = patient.RightSubmandibular
    lSubMasks = patient.LeftSubmandibularMasks
    rSubMasks = patient.RightSubmandibularMasks

    lSub_SUVs = np.zeros((9,2))
    rSub_SUVs = np.zeros((9,2)) #second dimension to hold number of points used for average, which will be removed after being used at the end
    print("Calculating Submandibular SUVs.")


    for slice_idx in range(0,pet_array.shape[1]):
        print("On Slice #" + str(slice_idx+1) + " for whole submandibular glands.")
        if np.amax(lSubMasks[0,slice_idx,:,:]) == 0 and np.amax(rSubMasks[0,slice_idx,:,:]) == 0: #if no contour on slice skip
            continue
        for x_idx in range(512):
            for y_idx in range(512):
                if lSubMasks[0,slice_idx, x_idx,y_idx] == 1:
                    lSub_SUVs[0,0] = lSub_SUVs[0,0] +  pet_array[0,slice_idx,x_idx,y_idx]
                    lSub_SUVs[0,1] = lSub_SUVs[0,1] + 1

                if rSubMasks[0,slice_idx, x_idx,y_idx] == 1:
                    rSub_SUVs[0,0] = rSub_SUVs[0,0] + pet_array[0,slice_idx,x_idx,y_idx]
                    rSub_SUVs[0,1] = rSub_SUVs[0,1] + 1    

        
    lSub_SUVs[0,0] = lSub_SUVs[0,0] / lSub_SUVs[0,1]
    rSub_SUVs[0,0] = rSub_SUVs[0,0] / rSub_SUVs[0,1]



    #Now calculate for the subsegments. 
    lSub8 = lSub.segmentedContours8
    rSub8 = rSub.segmentedContours8

    #check if subsegment directory exists and if not create one
    subsegDir = os.path.join(patient.path, "SubSegMasks")
    if not os.path.isdir(subsegDir):
        os.mkdir(subsegDir)

    print("Finished calculating SUVs for whole glands.")
    print("Beginning subsegment SUV calculations.")

    #first left par subsegs
    for subseg_idx, subsegment in enumerate(lSub8):
    
        subSegMasks_l = GetImageData.GetContourMasks(subsegment, pet_array)
        print("Calculating average SUV in subsegment " + str(int(subseg_idx+1))+ " of left submandibular.")
        for slice_idx in range(0,pet_array.shape[1]):
            
            if np.amax(subSegMasks_l[0,slice_idx,:,:]) == 0: #if no contour on slice skip
                continue

            layerBoolMask = subSegMasks_l[0,slice_idx,:,:] > 0
            subsegSUVs = pet_array[0,slice_idx,:,:][layerBoolMask]
            suv_sum = np.sum(subsegSUVs)
            numPoints = len(subsegSUVs)
            lSub_SUVs[subseg_idx+1, 0] = lSub_SUVs[subseg_idx+1, 0] + suv_sum
            lSub_SUVs[subseg_idx+1,1] = lSub_SUVs[subseg_idx+1,1] + numPoints

            
        lSub_SUVs[subseg_idx+1,0] = lSub_SUVs[subseg_idx+1,0] / lSub_SUVs[subseg_idx+1,1]
    lSub_SUVs = lSub_SUVs[:,0]

    #Now right par subsegs
    for subseg_idx, subsegment in enumerate(rSub8):

        subSegMasks_r = GetImageData.GetContourMasks(subsegment, pet_array)
        print("Calculating average SUV in subsegment " + str(int(subseg_idx+1))+ " of right submandibular.")

        for slice_idx in range(0,pet_array.shape[1]):
            if np.amax(subSegMasks_r[0,slice_idx,:,:]) == 0: #if no contour on slice skip
                continue

            layerBoolMask = subSegMasks_r[0,slice_idx,:,:] > 0
            subsegSUVs = pet_array[0,slice_idx,:,:][layerBoolMask]
            suv_sum = np.sum(subsegSUVs)
            numPoints = len(subsegSUVs)
            rSub_SUVs[subseg_idx+1, 0] = rSub_SUVs[subseg_idx+1, 0] + suv_sum
            rSub_SUVs[subseg_idx+1,1] = rSub_SUVs[subseg_idx+1,1] + numPoints

        

        rSub_SUVs[subseg_idx+1,0] = rSub_SUVs[subseg_idx+1,0] / rSub_SUVs[subseg_idx+1,1]
    rSub_SUVs = rSub_SUVs[:,0]


    rightSpearman = SpearmansRankCorrelation(rSub_SUVs[1:].tolist())
    leftSpearman = SpearmansRankCorrelation(lSub_SUVs[1:].tolist())
    rightSpearman_t = Get_t_value(rightSpearman, 8)
    leftSpearman_t = Get_t_value(leftSpearman, 8)
    rightSpearman_p = scipy.stats.t.sf(np.abs(rightSpearman_t), 6)
    leftSpearman_p = scipy.stats.t.sf(np.abs(leftSpearman_t), 6)
    rightPearson = PearsonsRankCorrelation(rSub_SUVs[1:].tolist())
    leftPearson = PearsonsRankCorrelation(lSub_SUVs[1:].tolist())
    rightPearson_t = Get_t_value(rightPearson, 8)
    leftPearson_t = Get_t_value(leftPearson, 8)
    rightPearson_p = scipy.stats.t.sf(np.abs(rightPearson_t), 6)
    leftPearson_p = scipy.stats.t.sf(np.abs(leftPearson_t), 6)
    #Now save these SUV stats into a csv

    csvPath = os.path.join(patient.path, "submandibular_stats.csv")

    with open(csvPath, 'w') as csvFile:
        filewriter = csv.writer(csvFile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['-', 'Left-Parotid', 'Right-Parotid'])
        filewriter.writerow(["whole-gland", lSub_SUVs[0], rSub_SUVs[0]]) 
        for idx in range(1, len(lSub_SUVs)):
            filewriter.writerow([str(int(idx)), lSub_SUVs[idx], rSub_SUVs[idx]]) 
        filewriter.writerow(["Spearmans Coefficient", leftSpearman, rightSpearman])
        filewriter.writerow(["Pearsons Coefficient", leftPearson, rightPearson])
        filewriter.writerow(["Spearman Significance", leftSpearman_p, rightSpearman_p])
        filewriter.writerow(["Pearson Significance", leftPearson_p, rightPearson_p])

    return [lSub_SUVs, rSub_SUVs]     

def ParotidSUVAnalysis(patient : Patient, stat="mean"):
    #this function takes a pet suv array and a structure and computes the mean suv within the structure    
    
    pet_array = patient.PETArray 
    lPar = patient.LeftParotid
    rPar = patient.RightParotid
    lParMasks = patient.LeftParotidMasks
    rParMasks = patient.RightParotidMasks

    lPar_SUVs = np.zeros((19,2))
    rPar_SUVs = np.zeros((19,2)) #second dimension to hold number of points used for average, which will be removed after being used at the end

    #also have a list for holding all voxel values to calculate percentiles with
    lPar_SUV_vals = []
    rPar_SUV_vals = []
    print("Calculating Parotid SUVs.")

    for slice_idx in range(0,pet_array.shape[1]):
        print("On Slice #" + str(slice_idx+1) + " for whole parotid glands.")
        if np.amax(lParMasks[0,slice_idx,:,:]) == 0 and np.amax(rParMasks[0,slice_idx,:,:]) == 0: #if no contour on slice skip
            continue
        for x_idx in range(512):
            for y_idx in range(512):
                if lParMasks[0,slice_idx, x_idx,y_idx] == 1:
                    lPar_SUVs[0,0] = lPar_SUVs[0,0] +  pet_array[0,slice_idx,x_idx,y_idx]
                    lPar_SUVs[0,1] = lPar_SUVs[0,1] + 1
                    lPar_SUV_vals.append(pet_array[0,slice_idx,x_idx,y_idx])

                if rParMasks[0,slice_idx, x_idx,y_idx] == 1:
                    rPar_SUVs[0,0] = rPar_SUVs[0,0] + pet_array[0,slice_idx,x_idx,y_idx]
                    rPar_SUVs[0,1] = rPar_SUVs[0,1] + 1    
                    rPar_SUV_vals.append(pet_array[0,slice_idx,x_idx,y_idx])

    if stat=="mean":    
        lPar_SUVs[0,0] = lPar_SUVs[0,0] / lPar_SUVs[0,1]
        rPar_SUVs[0,0] = rPar_SUVs[0,0] / rPar_SUVs[0,1]
    elif stat=="95":
        lPar_SUVs[0,0] = np.percentile(np.array(lPar_SUV_vals), 95)
        rPar_SUVs[0,0] = np.percentile(np.array(rPar_SUV_vals), 95)
    elif stat=="5":  
        lPar_SUVs[0,0] = np.percentile(np.array(lPar_SUV_vals), 5)
        rPar_SUVs[0,0] = np.percentile(np.array(rPar_SUV_vals), 5)  


    #Now calculate for the subsegments. 
    lPar18 = patient.LeftParotid.segmentedContours18
    rPar18 = patient.RightParotid.segmentedContours18

    #check if subsegment directory exists and if not create one
    subsegDir = os.path.join(patient.path, "SubSegMasks")
    if not os.path.isdir(subsegDir):
        os.mkdir(subsegDir)

    print("Finished calculating SUVs for whole glands.")
    print("Beginning subsegment SUV calculations.")

    #first left par subsegs
    for subseg_idx, subsegment in enumerate(lPar18):
        try: #try to load the masks for the subsegment

            with open(os.path.join(subsegDir, str(patient.name + "_left_subsegMasks_" + str(int(subseg_idx+1)) + ".txt")), "rb") as fp:
                subSegMasks_l = pickle.load(fp)
        except:      
            subSegMasks_l = GetImageData.GetContourMasks(subsegment, pet_array)
            # with open(os.path.join(subsegDir, str(patient.name + "_left_subsegMasks_" + str(int(subseg_idx+1)) + ".txt")), "wb") as fp:
            #     pickle.dump(subSegMasks_l, fp)
            ##takes too long to save
            lPar_SUV_vals = []
        print("Calculating average SUV in subsegment " + str(int(subseg_idx+1))+ " of left parotid.")
        for slice_idx in range(0,pet_array.shape[1]):
            
            if np.amax(subSegMasks_l[0,slice_idx,:,:]) == 0: #if no contour on slice skip
                continue

            layerBoolMask = subSegMasks_l[0,slice_idx,:,:] > 0
            subsegSUVs = pet_array[0,slice_idx,:,:][layerBoolMask]
            suv_sum = np.sum(subsegSUVs)

            lPar_SUV_vals.extend(subsegSUVs)

            numPoints = len(subsegSUVs)
            lPar_SUVs[subseg_idx+1, 0] = lPar_SUVs[subseg_idx+1, 0] + suv_sum
            lPar_SUVs[subseg_idx+1,1] = lPar_SUVs[subseg_idx+1,1] + numPoints

        if stat=="mean":    
            lPar_SUVs[subseg_idx+1,0] = lPar_SUVs[subseg_idx+1,0] / lPar_SUVs[subseg_idx+1,1]
        elif stat=="95":
            lPar_SUVs[subseg_idx+1,0] = np.percentile(np.array(lPar_SUV_vals), 95)
        elif stat=="5":  
            lPar_SUVs[subseg_idx+1,0] = np.percentile(np.array(lPar_SUV_vals), 5)
    lPar_SUVs = lPar_SUVs[:,0]

    #Now right par subsegs
    for subseg_idx, subsegment in enumerate(rPar18):
        try: #try to load the masks for the subsegment

            with open(os.path.join(subsegDir, str(patient.name + "_right_subsegMasks_" + str(int(subseg_idx+1)) + ".txt")), "rb") as fp:
                subSegMasks_r = pickle.load(fp)
        except:      
            subSegMasks_r = GetImageData.GetContourMasks(subsegment, pet_array)
            # with open(os.path.join(subsegDir, str(patient.name + "_right_subsegMasks_" + str(int(subseg_idx+1)) + ".txt")), "wb") as fp:
            #     pickle.dump(subSegMasks_r, fp)
            rPar_SUV_vals = []
        print("Calculating average SUV in subsegment " + str(int(subseg_idx+1))+ " of right parotid.")
        for slice_idx in range(0,pet_array.shape[1]):
            if np.amax(subSegMasks_r[0,slice_idx,:,:]) == 0: #if no contour on slice skip
                continue

            layerBoolMask = subSegMasks_r[0,slice_idx,:,:] > 0
            subsegSUVs = pet_array[0,slice_idx,:,:][layerBoolMask]
            suv_sum = np.sum(subsegSUVs)
            rPar_SUV_vals.extend(subsegSUVs)
            numPoints = len(subsegSUVs)
            rPar_SUVs[subseg_idx+1, 0] = rPar_SUVs[subseg_idx+1, 0] + suv_sum
            rPar_SUVs[subseg_idx+1,1] = rPar_SUVs[subseg_idx+1,1] + numPoints

        if stat=="mean":    
            rPar_SUVs[subseg_idx+1,0] = rPar_SUVs[subseg_idx+1,0] / rPar_SUVs[subseg_idx+1,1]
        elif stat=="95":
            rPar_SUVs[subseg_idx+1,0] = np.percentile(np.array(rPar_SUV_vals), 95)
        elif stat=="5":  
            rPar_SUVs[subseg_idx+1,0] = np.percentile(np.array(rPar_SUV_vals), 5)

    rPar_SUVs = rPar_SUVs[:,0]


    rightSpearman = SpearmansRankCorrelation(rPar_SUVs[1:].tolist())
    leftSpearman = SpearmansRankCorrelation(lPar_SUVs[1:].tolist())
    rightSpearman_t = Get_t_value(rightSpearman, 18)
    leftSpearman_t = Get_t_value(leftSpearman, 18)
    rightSpearman_p = scipy.stats.t.sf(np.abs(rightSpearman_t), 16)
    leftSpearman_p = scipy.stats.t.sf(np.abs(leftSpearman_t), 16)
    rightPearson = PearsonsRankCorrelation(rPar_SUVs[1:].tolist())
    leftPearson = PearsonsRankCorrelation(lPar_SUVs[1:].tolist())
    rightPearson_t = Get_t_value(rightPearson, 18)
    leftPearson_t = Get_t_value(leftPearson, 18)
    rightPearson_p = scipy.stats.t.sf(np.abs(rightPearson_t), 16)
    leftPearson_p = scipy.stats.t.sf(np.abs(leftPearson_t), 16)
    #Now save these SUV stats into a csv
    if stat=="mean":    
        csvPath = os.path.join(patient.path, "parotid_means.csv")
    elif stat=="95":
        csvPath = os.path.join(patient.path, "parotid_95_percentile.csv")
    elif stat=="5":  
        csvPath = os.path.join(patient.path, "parotid_5_percentile.csv")
    

    with open(csvPath, 'w') as csvFile:
        filewriter = csv.writer(csvFile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['-', 'Left-Parotid', 'Right-Parotid'])
        filewriter.writerow(["whole-gland", lPar_SUVs[0], rPar_SUVs[0]]) 
        for idx in range(1, len(lPar_SUVs)):
            filewriter.writerow([str(int(idx)), lPar_SUVs[idx], rPar_SUVs[idx]]) 
        filewriter.writerow(["Spearmans Coefficient", leftSpearman, rightSpearman])
        filewriter.writerow(["Pearsons Coefficient", leftPearson, rightPearson])
        filewriter.writerow(["Spearman Significance", leftSpearman_p, rightSpearman_p])
        filewriter.writerow(["Pearson Significance", leftPearson_p, rightPearson_p])

    return [lPar_SUVs, rPar_SUVs]    

def Population_Metrics_Parotid(stat="mean"):
    #this function computes the average, max and min suv inside each subsegment and whole gland of the left and parotid gland over all patients

    suvAvgs = [] #left is left, right is right parotid. thee first index is whole gland, followed by each subsegment
    suvMins = []
    suvMaxs = []
    suvSTDs = []
    suv95s = []
    suv5s = []

    for i in range(19):
        suvAvgs.append([0,0])
        suv95s.append([0,0])
        suv5s.append([0,0])
        suvMins.append([1e6, 1e6])
        suvMaxs.append([1e-6,1e-6])
        suvSTDs.append([0,0])

    #loop through first and get average, min and max suvs
    for i in range(1,31):
        patient_idx = "{:02d}".format(i)
#    print("Loading Patient: " + patient_idx)
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , patient_idx)
        if stat=="mean":
            csvPath = os.path.join(patientPath, 'parotid_stats.csv')
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
        elif stat=="95":
            csvPath = os.path.join(patientPath, 'parotid_95_percentile.csv')
            with open(csvPath) as csvFile:   
                csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
                for row in range(1, 20):          
                    for idx, rowData in csvData:

                        if int(idx) == int(row):
                            suv95s[row-1][0] = suv95s[row-1][0] + float(rowData[1])
                            suv95s[row-1][1] = suv95s[row-1][1] + float(rowData[2])
                            if float(rowData[1]) < suvMins[row-1][0]:
                                suvMins[row-1][0] = float(rowData[1])
                            if float(rowData[2]) < suvMins[row-1][1]:
                                suvMins[row-1][1] = float(rowData[2])   
        elif stat=="5":
            csvPath = os.path.join(patientPath, 'parotid_5_percentile.csv')
            with open(csvPath) as csvFile:   
                csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
                for row in range(1, 20):          
                    for idx, rowData in csvData:

                        if int(idx) == int(row):
                            suv95s[row-1][0] = suv95s[row-1][0] + float(rowData[1])
                            suv95s[row-1][1] = suv95s[row-1][1] + float(rowData[2])
                            if float(rowData[1]) < suvMins[row-1][0]:
                                suvMins[row-1][0] = float(rowData[1])
                            if float(rowData[2]) < suvMins[row-1][1]:
                                suvMins[row-1][1] = float(rowData[2])                                          
    for row in range(1, 20):
        suvAvgs[row-1][0] = suvAvgs[row-1][0] / 30
        suvAvgs[row-1][1] = suvAvgs[row-1][1] / 30
        suv95s[row-1][0] = suv95s[row-1][0] / 30
        suv95s[row-1][1] = suv95s[row-1][1] / 30
        suv5s[row-1][0] = suv5s[row-1][0] / 30
        suv5s[row-1][1] = suv5s[row-1][1] / 30


    #now loop back through and get stds.    
    for i in range(1,31):
        patient_idx = "{:02d}".format(i)
#    print("Loading Patient: " + patient_idx)
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , patient_idx)
        if stat=="mean":
            csvPath = os.path.join(patientPath, 'parotid_stats.csv')
        elif stat=="95":
            csvPath = os.path.join(patientPath, 'parotid_95_percentile.csv')
        elif stat=="5":
            csvPath = os.path.join(patientPath, 'parotid_5_percentile.csv')        

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
        patient_idx = "{:02d}".format(i)
#    print("Loading Patient: " + patient_idx)
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , patient_idx)
        if stat=="mean":
            csvPath = os.path.join(patientPath, 'parotid_stats.csv')
        elif stat=="95":
            csvPath = os.path.join(patientPath, 'parotid_95_percentile.csv')
        elif stat=="5":
            csvPath = os.path.join(patientPath, 'parotid_5_percentile.csv')        
        
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
    if stat=="mean":
        leftSubsegs, rightSubsegs = map(list, zip(*suvAvgs))
    elif stat=="95":
        leftSubsegs, rightSubsegs = map(list, zip(*suv95s))
    elif stat=="5":
        leftSubsegs, rightSubsegs = map(list, zip(*suv5s))
    leftSubsegs = leftSubsegs[1:]
    rightSubsegs = rightSubsegs[1:]
    rightSpearman = SpearmansRankCorrelation(CloneList(rightSubsegs))
    leftSpearman = SpearmansRankCorrelation(CloneList(leftSubsegs))
    rightSpearman_t = Get_t_value(rightSpearman, 18)
    leftSpearman_t = Get_t_value(leftSpearman, 18)
    rightSpearman_p = scipy.stats.t.sf(np.abs(rightSpearman_t), 16)
    leftSpearman_p = scipy.stats.t.sf(np.abs(leftSpearman_t), 16)
    rightPearson = PearsonsRankCorrelation(rightSubsegs)
    leftPearson = PearsonsRankCorrelation(leftSubsegs)
    rightPearson_t = Get_t_value(rightPearson, 18)
    leftPearson_t = Get_t_value(leftPearson, 18)
    rightPearson_p = scipy.stats.t.sf(np.abs(rightPearson_t), 16)
    leftPearson_p = scipy.stats.t.sf(np.abs(leftPearson_t), 16)
    print("Calculated SUV Stats")  
    
    #Now need to save these stats into the suv_stats csv file.
    path = os.path.join(parentDirectory, str("Statistics/PopulationParotidStats_" + stat + ".csv"))
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
def Population_Metrics_SM():
    #this function computes the average, max and min suv inside each subsegment and whole gland of the left and parotid gland over all patients

    suvAvgs = [] #left is left, right is right parotid. thee first index is whole gland, followed by each subsegment
    suvMins = []
    suvMaxs = []
    suvSTDs = []

    for i in range(9):
        suvAvgs.append([0,0])
        suvMins.append([1e6, 1e6])
        suvMaxs.append([1e-6,1e-6])
        suvSTDs.append([0,0])

    #loop through first and get average, min and max suvs
    for i in range(1,31):
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , str(i))
        csvPath = os.path.join(patientPath, 'submandibular_stats.csv')
        with open(csvPath) as csvFile:   
            csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
            for row in range(1, 10):          
                for idx, rowData in csvData:

                    if int(idx) == int(row):
                        suvAvgs[row-1][0] = suvAvgs[row-1][0] + float(rowData[1])
                        suvAvgs[row-1][1] = suvAvgs[row-1][1] + float(rowData[2])
                        if float(rowData[1]) < suvMins[row-1][0]:
                            suvMins[row-1][0] = float(rowData[1])
                        if float(rowData[2]) < suvMins[row-1][1]:
                            suvMins[row-1][1] = float(rowData[2])    
    for row in range(1, 10):
        suvAvgs[row-1][0] = suvAvgs[row-1][0] / 30
        suvAvgs[row-1][1] = suvAvgs[row-1][1] / 30


    #now loop back through and get stds.    
    for i in range(1,31):
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , str(i))
        csvPath = os.path.join(patientPath, 'submandibular_stats.csv')
        with open(csvPath) as csvFile:   
            csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
            for row in range(1, 10):          
                for idx, rowData in csvData:

                    if int(idx) == int(row):
                        suvSTDs[row-1][0] = suvSTDs[row-1][0] + (float(rowData[1])-suvAvgs[row-1][0])**2
                        suvSTDs[row-1][1] = suvSTDs[row-1][1] + (float(rowData[2])-suvAvgs[row-1][1])**2
    for row in range(1, 10):
        suvSTDs[row-1][0] = (suvSTDs[row-1][0] / 30)**0.5
        suvSTDs[row-1][1] = (suvSTDs[row-1][1] / 30)**0.5
    


    #Now get the number of significant correlations (spearman and pearson) using a significant cutoff of p = 0.01, 0.02, 0.05, 0.1 for left and right parotids
    numSig_Spear_Left = [0,0,0,0]
    numSig_Spear_Right = [0,0,0,0]
    numSig_Pears_Left = [0,0,0,0]
    numSig_Pears_Right = [0,0,0,0]
    for i in range(1,31):
        patientPath = os.path.join(parentDirectory, "SG_PETRT" , str(i))
        csvPath = os.path.join(patientPath, 'submandibular_stats.csv')
        
        with open(csvPath) as csvFile:   
            csvData = list(enumerate(csv.reader(csvFile, delimiter=',', quotechar='|')))       
            leftSpear, rightSpear = [float(i) for i in csvData[12][1][1:]]
            leftPears, rightPears = [float(i) for i in csvData[13][1][1:]]
            
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
    rightSpearman_t = Get_t_value(rightSpearman, 8)
    leftSpearman_t = Get_t_value(leftSpearman, 8)
    rightSpearman_p = scipy.stats.t.sf(np.abs(rightSpearman_t), 6)
    leftSpearman_p = scipy.stats.t.sf(np.abs(leftSpearman_t), 6)
    rightPearson = PearsonsRankCorrelation(rightSubseg_avgs)
    leftPearson = PearsonsRankCorrelation(leftSubseg_avgs)
    rightPearson_t = Get_t_value(rightPearson, 8)
    leftPearson_t = Get_t_value(leftPearson, 8)
    rightPearson_p = scipy.stats.t.sf(np.abs(rightPearson_t), 6)
    leftPearson_p = scipy.stats.t.sf(np.abs(leftPearson_t), 6)
    print("Calculated SUV Stats")  
    
    #Now need to save these stats into the suv_stats csv file.
    path = os.path.join(parentDirectory, "Statistics/PopulationSubmandibularStats.csv")
    with open(path, 'w') as csvFile:
        filewriter = csv.writer(csvFile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['Subsegment', 'Left-Sub','Left-Sub-SD', 'Right-Sub', 'Right-Sub-SD'])
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
    elif len(suvs) == 8:    
        importanceVals = SubImportance()
    importanceRanks = GetListRank(importanceVals)      
    suvRanks = GetListRank(suvs)
    n = len(suvRanks)
    sumRankDifs = 0
    for i in range(len(suvRanks)):
        sumRankDifs = sumRankDifs + (importanceRanks[i] - suvRanks[i])**2
    rho = 1 - ((6*sumRankDifs)/(n*(n**2-1)))    
    return rho
def PearsonsRankCorrelation(suvs, gland="parotid"):
    #This function computes the spearmans rank coefficient for a list of suv values

    if len(suvs) == 18: #if 18 subsegments
        importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
                1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
                0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
                0.272984146341463,   0.0484897103658537,  0.035493902439024]
    elif len(suvs) == 8:    
        importanceVals = SubImportance()

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

    SUV_Metrics()                
