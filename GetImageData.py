from math import nan
import os
from numpy.lib.function_base import average 
import pydicom
import numpy as np
from fastDamerauLevenshtein import damerauLevenshtein
from shapely.geometry.polygon import Polygon
from Contours import Contours
import Chopper
from operator import delitem, itemgetter
from PIL import Image, ImageDraw
from Patient import Patient
import pickle 
import statistics
import copy
import csv
import scipy.stats
from Data_Analyzing import Get_t_value, SpearmansRankCorrelation, GetListRank, PearsonsRankCorrelation, CloneList, DicomSaver



parentDirectory = os.getcwd()

def GetAvgAge():
    path = os.path.join(parentDirectory, "SG_PETRT")
    patients = os.listdir(path)
    ages = []
    for patient in patients:
        patient_path = os.path.join(path, patient)
        folders = os.listdir(patient_path)
        for folder in folders:
            if "PET_" in folder:
                pet_path = os.path.join(patient_path, folder)
                pet_files = os.listdir(pet_path)
                pet_file = os.path.join(pet_path, pet_files[0])
                metadata = pydicom.dcmread(pet_file)
                age = metadata[0x0010, 0x1010].value
                ages.append(int(age[0:3]))

                break
    avg_age = statistics.mean(ages)
    min_age = min([i for i in ages if i != 0])
    max_age = max(ages)
    deviation = max(avg_age-min_age, max_age-avg_age)
    print("Average age for the cohort of " + str(len(ages)) + " patients is " + str(avg_age) + " +/- " + str(deviation)   )  

def GetAvgWeight():
    path = os.path.join(parentDirectory, "SG_PETRT")
    patients = os.listdir(path)
    weights = []
    for patient in patients:
        patient_path = os.path.join(path, patient)
        folders = os.listdir(patient_path)
        for folder in folders:
            if "PET_" in folder:
                pet_path = os.path.join(patient_path, folder)
                pet_files = os.listdir(pet_path)
                pet_file = os.path.join(pet_path, pet_files[0])
                metadata = pydicom.dcmread(pet_file)
                weight = metadata[0x0010, 0x1030].value
                weights.append(int(weight))


                break
    avg_weight = statistics.mean(weights)
    min_weight = min(weights)
    max_weight = max(weights)
    deviation = max(avg_weight-min_weight, max_weight-avg_weight)
    print("Average weight for the cohort of " + str(len(weights)) + " patients is " + str(avg_weight) + " +/- " + str(deviation)  )  

def GetSexes():
    path = os.path.join(parentDirectory, "SG_PETRT")
    patients = os.listdir(path)
    M = 0
    F = 0
    for patient in patients:
        patient_path = os.path.join(path, patient)
        folders = os.listdir(patient_path)
        for folder in folders:
            if "PET_" in folder:
                pet_path = os.path.join(patient_path, folder)
                pet_files = os.listdir(pet_path)
                pet_file = os.path.join(pet_path, pet_files[0])
                metadata = pydicom.dcmread(pet_file)
                sex = metadata[0x0010, 0x0040].value
                if sex == 'M':
                    M += 1
                elif sex == 'F':
                    F += 1
                else:
                    print("")    



                break
    print("sexes for the cohort: " + str(M) + "M / " + str(F) + "F"    )        



def RemoveIslands(contours):
    contours = sorted(contours, key=lambda x: x[0][2])
    indices_to_remove = []
    
    zVals = []
    for i in range(len(contours)):
        #first loop through and get distance between points statistics
        distances = []
        zVals.append(contours[i][0][2])
    
    duplicateSlices = GetDuplicateSlices(contours)

    for slice, indices in duplicateSlices.items():
        #first check to see if any are tiny contours:
        lengths = []
        for index in indices:
            lengths.append(len(contours[index]))
        longestLength = max(lengths)
        longestIndex = lengths.index(longestLength)
        #now see if this is more than 3x larger than all other contours, and if it is, return it
        lengths.pop(longestIndex)  
        muchLarger = True  
        for value in lengths:
            if longestLength / value < 3:
                muchLarger = False

        if muchLarger == True:
            indices.pop(longestIndex)
            for item in indices:
                indices_to_remove.append(item)
            continue
                 
        #first get average x and y of slices immediately above and below (or only one if at top or bottom)
        if min(indices) == 0:
            above_idx = max(indices)+1
            xVals = []
            yVals = []
            for point in contours[above_idx]:
                xVals.append(point[0])
                yVals.append(point[1])
            averageX = average(xVals)
            averageY = average(yVals)

        elif max(indices) == len(contours)-1:
            below_idx = min(indices)-1
            xVals = []
            yVals = []
            for point in contours[below_idx]:
                xVals.append(point[0])
                yVals.append(point[1])
            averageX = average(xVals)
            averageY = average(yVals)
        else:
            above_idx = max(indices)+1
            below_idx = min(indices)-1
            xVals = []
            yVals = []
            for point in contours[below_idx]:
                xVals.append(point[0])
                yVals.append(point[1])
            for point in contours[above_idx]:
                xVals.append(point[0])
                yVals.append(point[1])    
            averageX = average(xVals)
            averageY = average(yVals)

        #Now get the slice statistics: 
        displacements = [] #holds the displacement of each slice from averages below/on top , and we return the one closest to averages
        for index in indices:
            xVals = []
            yVals = []
            for point in contours[index]:
                xVals.append(point[0])
                yVals.append(point[1])
            averageX_slice = average(xVals)
            averageY_slice = average(yVals)

            displacement = ( (averageX-averageX_slice)**2 + (averageY-averageY_slice)**2)    
            displacements.append(displacement)
        minDisplacement = min(displacements)


        removeIndices = indices.copy()
        removeIndices.pop(displacements.index(minDisplacement)) #remove the one we want to keep from the remove list. All removeIndices will be taken out of contours. 
        for item in removeIndices:
            indices_to_remove.append(item)
    
    indices_to_remove.sort(reverse=True)
    for item in indices_to_remove:
        contours.pop(item)    
    return contours


  

def GetDuplicateSlices(contours):
    #this function takes a contour list and returns indices of contours which have same z value in a dictionary
    duplicateSliceDict = dict()
    for i, slice in enumerate(contours):
        if slice[0][2] in duplicateSliceDict:
            duplicateSliceDict[slice[0][2]].append(i)
        else:
            duplicateSliceDict[slice[0][2]] = [i]
    #filter out z values with only one slice
    duplicateSliceDict = {key:value for key, value in duplicateSliceDict.items() if len(value) > 1}            
    return duplicateSliceDict




def GetContours(patientPath, organ, subsegmentation = [2,1,2]): 
    patientItems = os.listdir(patientPath)
    structFiles = []
    for item in patientItems:
        if "STRUCT" in item:
            structFiles.append(os.path.join(patientPath, item))
    
    structure, structureROINum, struct_idx= FindStructure(structFiles, organ)  
    structsMeta = pydicom.dcmread(structFiles[struct_idx]).data_element("ROIContourSequence")
        

    if structureROINum == 1111:
        return nan
    contourList = []
    for contourInfo in structsMeta:
        if contourInfo.get("ReferencedROINumber") == structureROINum: #get the correct contour for the given organ
            for contoursequence in contourInfo.ContourSequence: 
                contourList.append(contoursequence.ContourData)
                #But this is easier to work with if we convert from a 1d to a 2d list for contours ( [ [x1,y1,z1], [x2,y2,z2] ... ] )
                numContourPoints = 0 
                tempContour = []
                i = 0
                while i < len(contourList[-1]):
                    x = float(contourList[-1][i])
                    y = float(contourList[-1][i + 1])
                    z = float(contourList[-1][i + 2])
                    tempContour.append([x, y, z ])
                    i += 3
                    if (numContourPoints > 0):
                        pointListy = np.vstack((pointListy, np.array([x,y,z])))
                    else:
                        pointListy = np.array([x,y,z])   
                    numContourPoints+=1         
                contourList[-1] = tempContour 
    contourList = RemoveIslands(CloneList(contourList))            
    contours = Contours(organ, structure, contourList)   
    
    Chopper.OrganChopper(contours, subsegmentation, organ) #Get 18ths subsegments
    
    return contours             


       
def GetPatientName(patient):
    patientFolders = os.listdir(patient)
    for folder in patientFolders:
        if "PET" in folder:
            petFolder = os.path.join(patient, folder)
            petImages = os.listdir(petFolder)

    imagesandSlices = []
    for idx, image in enumerate(petImages):
        imagePath = os.path.join(petFolder, image)
        metadata = pydicom.dcmread(imagePath)
        return metadata[0x0010,0x0020].value
     
def GetCTArray(patient):
    patientFolders = os.listdir(patient)
    for folder in patientFolders:
        if "CT" in folder and "STRUCT" not in folder:
            ctFolder = os.path.join(patient, folder)
            ctImages = os.listdir(ctFolder)

    imagesandSlices = []
    for idx, image in enumerate(ctImages):
        imagePath = os.path.join(ctFolder, image)
        metadata = pydicom.dcmread(imagePath)
        pixelData = metadata.pixel_array
        sliceLocation = metadata[0x0020, 0x1041].value
        imagesandSlices.append([pixelData, sliceLocation])

        if idx == 0:
            xLen, yLen = pixelData.shape
            patientName = metadata[0x0010,0x0020].value
            ipp = metadata[0x0020,0x0032].value #image position patient
            iop = metadata[0x0020, 0x0037].value
            pixelSpacing = metadata[0x0028,0x0030].value
            array = np.zeros((4, len(ctImages), xLen, yLen))


        

    #now sort list according to slice location
    imagesandSlices.sort(key= lambda x: x[1])
    print("Collected CT data for:") 
    for idx, item in enumerate(imagesandSlices):
        array[0,idx,:,:] = item[0]
        array[3,idx,:,:] = item[1]
        
        for x_idx in range(xLen):
            for y_idx in range(yLen):
                x = ipp[0] + x_idx*pixelSpacing[0]
                y = ipp[0] + y_idx*pixelSpacing[1]
                array[1,idx,x_idx,y_idx] = x
                array[2,idx,x_idx,y_idx] = y
        print("    " + str(idx+1) + "/" + str(len(imagesandSlices)) + " DICOM files")         
    return array
    print("Successfully retrieved CT Array for " + patientName)               

def HHMMSS_To_S(val):
    s = 0
    hours = float(val[0:2]) * 60 * 60
    mins = float(val[2:4]) * 60
    seconds = float(val[4:])
    s = hours + mins + seconds
    return s

def GetSUVArray(patient):
    #This returns a 4d array with the first dimension having 4 channels - one for pet image, and other 3 for x,y,z coordinates of pixel

    patientFolders = os.listdir(patient)
    for folder in patientFolders:
        if "PET" in folder:
            petFolder = os.path.join(patient, folder)
            petImages = os.listdir(petFolder)

    imagesandSlices = []
    for idx, image in enumerate(petImages):
        imagePath = os.path.join(petFolder, image)
        metadata = pydicom.dcmread(imagePath)
        pixelData = metadata.pixel_array
        sliceLocation = metadata[0x0020, 0x1041].value
        rescaleSlope = float(metadata[0x0028, 0x1053].value)
        patientWeight = float(metadata[0x0010, 0x1030].value)
        acquisitionTime = metadata[0x0008,0x0031].value
        acquisitionTime = HHMMSS_To_S(acquisitionTime)
        halfLife = float(metadata[0x0054, 0x0016][0][0x0018, 0x1075].value)
        radStartTime = metadata[0x0054, 0x0016][0][0x0018, 0x1072].value
        radStartTime = HHMMSS_To_S(radStartTime)
        totalDose = float(metadata[0x0054,0x0016][0][0x0018, 0x1074].value)
        
        imagesandSlices.append([pixelData, sliceLocation, rescaleSlope, acquisitionTime])
        if idx == 0:
            xLen, yLen = pixelData.shape
            patientName = metadata[0x0010,0x0020].value
            ipp = metadata[0x0020,0x0032].value #image position patient
            iop = metadata[0x0020, 0x0037].value
            pixelSpacing = metadata[0x0028,0x0030].value
            array = np.zeros((5, len(petImages), xLen, yLen)) #first index in SUVBW, then second is BqmL, then position arrays (xyz)




    #now sort list according to slice location
    imagesandSlices.sort(key= lambda x: x[1])
    print("Collected PET SUV data for:") 
    for idx, item in enumerate(imagesandSlices):
        totalDoseCorrected = float(totalDose * 2 **((radStartTime-item[3])/halfLife))
        rescaleSlope = item[2]
        
        array[0,idx,:,:] = item[0] * rescaleSlope * patientWeight * 1000 / totalDoseCorrected
        array[1,idx,:,:] = item[0] * rescaleSlope
        array[4,idx,:,:] = item[1]
        
        for x_idx in range(xLen):
            for y_idx in range(yLen):
                x = ipp[0] + x_idx*pixelSpacing[0]
                y = ipp[1] + y_idx*pixelSpacing[1]
                array[2,idx,x_idx,y_idx] = x
                array[3,idx,x_idx,y_idx] = y
        print("    " + str(idx+1) + "/" + str(len(imagesandSlices)) + " DICOM files")         
    return array
    print("Got SUV Array for " + patientName)

def FindStructure(structureList, organ, invalidStructures = []):
    """Finds the matching structure to a given organ in a patient's
       dicom file metadata. 
     
    Args:
        structureList (List): a list of paths to RTSTRUCT files in the patient folder
        organ (str): the organ to find the matching structure for
        invaidStructures (list): a list of structures that the matching 
            structure cannot be, defaults to an empty list

    Returns: 
        str, int: the matching structure's name in the metadata, the 
            matching structure's ROI number in the metadata. Returns "", 
            1111 if no matching structure is found in the metadata
        
    """

    #Here we take the string for the desired structure (organ) and find the matching structure for each patient. 
    #The algorithm is to first make sure that the organ has a substring matching the ROI with at least 3 characters,
    #then out of the remaining possiblities, find top 3 closest fits with damerau levenshtein algorithm, then check to make sure that they are allowed to match according to rules defined in AllowedToMatch(). There should then ideally
    # be only one possible match, but if there are two, take the first in the list.   

    #Get a list of all structures in structure set: 
    unfilteredStructures = []

    for fileNum, file in enumerate(structureList):
        
        roiSequence = pydicom.dcmread(file).data_element("StructureSetROISequence")
        for element in roiSequence:
            if element.get("ROIName").lower() not in invalidStructures:
                #first check if need to do special matching for limbus name:
                if element.get("ROIName").lower() == "sm_l" and organ  == "Left Submandibular":
                    roiNumber = element.get("ROINumber")
                    return element.get("ROIName").lower(), roiNumber, fileNum
                if element.get("ROIName").lower() == "sm_r" and organ  == "Right Submandibular":
                    roiNumber = element.get("ROINumber")
                    return element.get("ROIName").lower(), roiNumber, fileNum    
                unfilteredStructures.append(element.get("ROIName").lower())           
    #Now find which is the best fit.
    #First filter out any structures without at least a 3 character substring in common
    structures = []
    for structure in unfilteredStructures:
        valid = True
        if LongestSubstring(structure, organ) < 3:
            valid = False
        if not AllowedToMatch(organ, structure): 
            valid = False
        #Add to structures if valid
        if valid:
            structures.append(structure)
    #Now test string closeness to find
    closestStrings = [["",100],["",100],["",100]] #has to be in the top 3 closest strings to check next conditions
    for structure in structures:
        closeness = StringDistance(structure, organ)
        closestStrings.sort(key=itemgetter(1)) #Sort by the closeness value, and not the structure names
        for value in range(len(closestStrings)):
            if closeness < closestStrings[value][1]: #If closer than a value already in the top 3
                closestStrings[value] = [structure, closeness]
                break
    
    if len(closestStrings) == 0:
        return "", 1111, 0    
    #Now return the organ that is remaining and has closest string
    fileNum = -1
    for file in structureList:
        fileNum = fileNum + 1
        roiSequence = pydicom.dcmread(file).data_element("StructureSetROISequence")
        for element in roiSequence:
            if element.get("ROIName").lower() == closestStrings[0][0]:
                roiNumber = element.get("ROINumber")
    try:
        return closestStrings[0][0], roiNumber, fileNum
    except:
        return "", 1111, 0 #error code for unfound match.                                                                                                          

def AllowedToMatch(s1, s2):
    """Determines whether or not s1 and s2 are allowed to match 
       based on if they both contain the correct substrings.
     
    Args:
        s1 (str): first string to determine match
        s2 (str): second string to determine match

    Returns: 
        allowed (bool): True if the strings are allowed to match, 
            false otherwise
        
    """

    s1 = s1.lower()
    s2 = s2.lower()
    allowed = True
    keywords = []
    #You can't have only one organ with one of these keywords...
    keywords.append("prv")
    keywords.append("tub")
    keywords.append("brain")
    keywords.append("ptv")
    keywords.append("stem")
    keywords.append("node")
    keywords.append("cord")
    keywords.append("chi")
    keywords.append("opt")
    keywords.append("oral")
    keywords.append("nerv")
    keywords.append("par")
    keywords.append("globe")
    keywords.append("lip")
    keywords.append("cav")
    keywords.append("sub")
    keywords.append("test")
    keywords.append("fact")
    #keywords can't be in only one of two string names: 
    for keyword in keywords:
        num = 0
        if keyword in s1:
            num += 1
        if keyword in s2:
            num += 1
        if num == 1:
            allowed = False        

    #Cant have left and no l in other, or right and no r
    if "left" in s1:
        if "l" not in s2:
            allowed = False      
    if "left" in s2:
        if "l" not in s1:
            allowed = False    
    #its tricky matching up left and right organs sometimes with all the conventions used... this makes sure that both are left or both are right
    if ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3]) or ("_l" == s1[-2:]):
        if not (("lpar" in s2) or ("lsub" in s2) or ("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("-lt " in s2) or ("lt " == s2[0:3]) or ("_l" == s2[-2:])):   
            allowed = False  
    if (("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("-lt " in s2)or ("lt " == s2[0:3]) or ("_l" == s2[-2:])):  
        if not (("lpar" in s1) or ("lsub" in s1) or ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3]) or ("_l" == s1[-2:])):
            allowed = False        
    
    if ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)or ("right" in s1) or ("_r" == s1[-2:]):
        if not (("rpar" in s2) or ("rsub" in s2) or ("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("-rt" in s2) or ("_r" == s2[-2:])):   
            allowed = False
    if (("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("-rt" in s2) or ("_r" == s2[-2:])): 
        if not (("rpar" in s1) or ("rsub" in s1) or ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1) or ("_r" == s1[-2:])):
            allowed = False
    return allowed


def StringDistance(s1, s2):
    """returns the Damerau-Levenshtein distance between two strings

    Args:
        s1 (string): string one which is to be compared with string 2.
        s2 (string): string two which is to be compared with string 1.

    Returns:
        (int): the Damerau Levenshtein distance between s1 and s2, which indicates how different the two strings are in terms of the amount of deletion, insertion, substitution, and transposition operations required to equate the two.

    """
    return damerauLevenshtein(s1,s2,similarity=False)

def LongestSubstring(s1,s2):
    """Finds the length of the longest substring that is in 
       both s1 and s2.
     
    Args:
        s1 (str): the first string to find the longest substring in
        s2 (str): the second string to find the longest substring in

    Returns: 
        longest (int): the length of the longest substring that is in 
            s1 and s2
        
    """

    m = len(s1)
    n = len(s2)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if s1[i] == s2[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(s1[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(s1[i-c+1:i+1])
    return longest  

# def FindStructure(metadata, organ, invalidStructures = []):
#     #Here we take the string for the desired structure (organ) and find the matching structure for each patient. 
#     #The algorithm is to first make sure that the organ has a substring matching the ROI with at least 3 characters,
#     #then out of the remaining possiblities, find top 3 closest fits with damerau levenshtein algorithm, then check to make sure that they are allowed to match according to rules defined in AllowedToMatch(). There should then ideally
#     # be only one possible match, but if there are two, take the first in the list.   
#     #Get a list of all structures in structure set: 
#     unfilteredStructures = []
#     for element in metadata:
#         if element.get("ROIName").lower() not in invalidStructures:
#             unfilteredStructures.append(element.get("ROIName").lower())           
#     #Now find which is the best fit.
#     #First filter out any structures without at least a 3 character substring in common
#     structures = []
#     for structure in unfilteredStructures:
#         valid = True
#         if LongestSubstring(structure, organ) < 3:
#             valid = False
#         if not AllowedToMatch(organ, structure): 
#             valid = False
#         #Add to structures if valid
#         if valid:
#             structures.append(structure)
#     #Now test string closeness to find
#     closestStrings = [["No Match",100],["No Match",100],["No Match",100]] #has to be in the top 3 closest strings to check next conditions
#     for structure in structures:
#         closeness = StringDistance(structure, organ)
#         closestStrings.sort(key=itemgetter(1)) #Sort by the closeness value, and not the structure names
#         for value in range(len(closestStrings)):
#             if closeness < closestStrings[value][1]: #If closer than a value already in the top 3
#                 closestStrings[value] = [structure, closeness]
#                 break
    
#     if closestStrings[0][1] == 100:
#         return "No Match", 1111    
#     #Now return the organ that is remaining and has closest string
#     for element in metadata:
#         if element.get("ROIName").lower() == closestStrings[0][0]:
#             roiNumber = element.get("ROINumber")
#     try:
#         return closestStrings[0][0], roiNumber 
#     except:
#         return "No Match", 1111 #error code for unfound match.    

# def AllowedToMatch(s1, s2):
#     s1 = s1.lower()
#     s2 = s2.lower()
#     allowed = True
#     keywords = []
#     #You can't have only one organ with one of these keywords...
#     keywords.append("prv")
#     keywords.append("brain")
#     keywords.append("ptv")
#     keywords.append("stem")
#     keywords.append("node")
#     keywords.append("cord")
#     keywords.append("chi")
#     keywords.append("opt")
#     keywords.append("oral")
#     keywords.append("nerv")
#     keywords.append("par")
#     keywords.append("globe")
#     keywords.append("lip")
#     keywords.append("cav")
#     keywords.append("sub")
#     keywords.append("test")
#     keywords.append("fact")
#     #keywords can't be in only one of two string names: 
#     for keyword in keywords:
#         num = 0
#         if keyword in s1:
#             num += 1
#         if keyword in s2:
#             num += 1
#         if num == 1:
#             allowed = False        

#     #Cant have left and no l in other, or rightt and no r
#     if "left" in s1:
#         if "l" not in s2:
#             allowed = False      
#     if "left" in s2:
#         if "l" not in s1:
#             allowed = False    
#     #its tricky matching up left and right organs sometimes with all the conventions used... this makes sure that both are left or both are right
#     if ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3]):
#         if not (("lpar" in s2) or ("lsub" in s2) or ("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("_lt" in s2) or ("-lt " in s2) or ("lt " == s2[0:3])):   
#             allowed = False  
#     if (("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("-lt " in s2)or ("lt " == s2[0:3])):  
#         if not (("lpar" in s1) or ("lsub" in s1) or ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3])):
#             allowed = False        
    
#     if ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)or ("right" in s1):
#         if not (("rpar" in s2) or ("rsub" in s2) or ("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("_rt" in s2) or ("-rt" in s2) ):   
#             allowed = False
#     if (("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("-rt" in s2) ): 
#         if not (("rpar" in s1) or ("rsub" in s1) or ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)):
#             allowed = False
#     return allowed


# def StringDistance(s1, s2):
#     return damerauLevenshtein(s1,s2,similarity=False)

#     return d[lenstr1-1,lenstr2-1]
# def LongestSubstring(s1,s2):
#     m = len(s1)
#     n = len(s2)
#     counter = [[0]*(n+1) for x in range(m+1)]
#     longest = 0
#     lcs_set = set()
#     for i in range(m):
#         for j in range(n):
#             if s1[i] == s2[j]:
#                 c = counter[i][j] + 1
#                 counter[i+1][j+1] = c
#                 if c > longest:
#                     lcs_set = set()
#                     longest = c
#                     lcs_set.add(s1[i-c+1:i+1])
#                 elif c == longest:
#                     lcs_set.add(s1[i-c+1:i+1])
#     return longest  


def GetContourMasks(contours, Array):
    numImagesSlices, xLen, yLen = Array.shape[1:]
    contourMasks = np.zeros((2,numImagesSlices,xLen,yLen)) #2 channels, one for filled and one for unfilled

    contours = CartesianToPixelCoordinates(CloneList(contours), Array)
    for idx in range(numImagesSlices):#loop through all slices creating a mask for the contours
        for contour in contours:
            if len(contour) < 3:
                continue
            if abs(int(round(contour[0][2], 2)*100) - int(round(Array[4,idx,0,0], 2)*100)) < 2: #if contour is on the current slice
                contourMaskFilled = Image.new('L', (xLen, yLen), 0 )
                contourMaskUnfilled = Image.new('L', (xLen, yLen), 0 )
                contourPoints = []

                for point in contour:
                    contourPoints.append((int(point[0]), int(point[1])))
                contourPolygon = Polygon(contourPoints)
                ImageDraw.Draw(contourMaskFilled).polygon(contourPoints, outline= 1, fill = 1)
                ImageDraw.Draw(contourMaskUnfilled).polygon(contourPoints, outline= 1, fill = 0)
                contourMaskFilled = np.array(contourMaskFilled)
                contourMaskUnfilled = np.array(contourMaskUnfilled)
                contourMasks[0, idx, :,:] = contourMaskFilled
                contourMasks[1,idx,:,:] = contourMaskUnfilled    
                break                     
    return contourMasks                


def CartesianToPixelCoordinates(contours, array):
    #convert x and y values for a contour into the pixel indices where they are on the pet array
    xVals = array[2,0,:,0]
    yVals = array[3,0,0,:]
    for contour in contours: 
        for point in contour:
            point[0] = min(range(len(xVals)), key=lambda i: abs(xVals[i]-point[0]))
            point[1] = min(range(len(yVals)), key=lambda i: abs(yVals[i]-point[1]))
    return contours        

def ImageUpsizer(array, newDimensions):
    """Supersizes an array by the factor given.   

    Args:
        array (2D numpy array): an array to be supersized
        newDimensions (list): the size of the array which is to be returned

    Returns:
        newArray (2D numpy array): the supersized array

    """
    xLen, yLen = array.shape
    xLen_new, yLen_new = newDimensions

    intermediateArray= np.zeros((xLen, yLen_new)) #first make an intermediate array with the sampling along y
    newArray = np.zeros((xLen_new, yLen_new))

    scaleSize = float(xLen_new / xLen)
    for i in range(xLen):
        for j in range(yLen_new):
            pixelFactor = float(j/scaleSize)
            leftPixel = int(pixelFactor)
            rightPixel = leftPixel + 1
            pixelFactor = pixelFactor - leftPixel
            if rightPixel < yLen:
                newVal = (1-pixelFactor) * array[i, leftPixel] + pixelFactor*array[i, rightPixel]
            elif rightPixel == yLen:
                newVal = array[i, leftPixel]    
            intermediateArray[i, j] = newVal

    for j in range(yLen_new):
        for i in range(xLen_new):
            pixelFactor = float(i/scaleSize)
            leftPixel = int(pixelFactor)
            rightPixel = leftPixel + 1
            pixelFactor = pixelFactor - leftPixel
            if rightPixel < xLen:
                newVal = (1-pixelFactor) * intermediateArray[leftPixel, j] + pixelFactor*intermediateArray[rightPixel, j]
            elif rightPixel == xLen:
                newVal = intermediateArray[leftPixel,j]
            newArray[i,j] = newVal

    return newArray


    # #Take an array and supersize it by the factor given
    # xLen, yLen = array.shape
    # newArray = np.zeros((factor * xLen, factor * yLen))
    # #first get the original values in to the grid: 
    # for i in range(xLen):
    #     for j in range(yLen):
    #         newArray[i * factor, j * factor] = array[i,j]
    # #sample along first dim
    # for j in range(yLen):
    #     for i in range(xLen - 1):
    #         insert = 1 
    #         while insert <= factor - 1:
    #             newArray[i * factor + insert, j * factor] = newArray[i * factor, j * factor] + (insert / factor) * (newArray[(i+1) * factor, j * factor]- newArray[i * factor, j * factor])
    #             insert += 1
    # #sample along second dim
    # for i in range(xLen * factor):
    #     for j in range(yLen - 1):
    #         insert = 1 
    #         while insert <= factor - 1:
    #             newArray[i, j * factor + insert] = newArray[i, j * factor] + (insert / factor) * (newArray[i, (j+1) * factor]- newArray[i, j * factor])
    #             insert += 1
    # return newArray









if __name__ == '__main__':
    #GetContours("/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1")
    #GetSUVArray("/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1")
    # SpearmansRankCorrelation([0.751310670731707,  0.526618902439024,   0.386310975609756,
    #         1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
    #         0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
    #         0.272984146341463,   0.0484897103658537,  0.035493902439024]) 
    GetAvgAge()