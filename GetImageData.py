from math import nan
import os 
import pydicom
import numpy as np
from fastDamerauLevenshtein import damerauLevenshtein
from shapely.geometry.polygon import Polygon
from Contours import Contours
import Chopper
from operator import itemgetter
from PIL import Image, ImageDraw


def GetContours(patientPath, organ): 
    patientItems = os.listdir(patientPath)
    for item in patientItems:
        if "STRUCT" in item:
            contourPath = os.path.join(patientPath, item)
    structsMeta = pydicom.dcmread(contourPath).data_element("ROIContourSequence")
    structure, structureROINum= FindStructure(pydicom.dcmread(contourPath).data_element("StructureSetROISequence"), organ)  
        

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
    contours = Contours(organ, structure, contourList)   
    Chopper.OrganChopper(contours, [2, 1, 2], organ) #Get 18ths subsegments
    
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
    for idx, item in enumerate(imagesandSlices):
        array[0,idx,:,:] = item[0]
        array[3,idx,:,:] = item[1]
        
        for x_idx in range(xLen):
            for y_idx in range(yLen):
                x = ipp[0] + x_idx*pixelSpacing[0]
                y = ipp[0] + y_idx*pixelSpacing[1]
                array[1,idx,x_idx,y_idx] = x
                array[2,idx,x_idx,y_idx] = y
        print("Collected CT data for " + str(idx+1) + "/" + str(len(imagesandSlices)) + " DICOM files")         
    return array
    print("Got CT Array for " + patientName)               


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
        imagesandSlices.append([pixelData, sliceLocation])

        if idx == 0:
            xLen, yLen = pixelData.shape
            patientName = metadata[0x0010,0x0020].value
            ipp = metadata[0x0020,0x0032].value #image position patient
            iop = metadata[0x0020, 0x0037].value
            pixelSpacing = metadata[0x0028,0x0030].value
            array = np.zeros((4, len(petImages), xLen, yLen))

            #Now get the x,y position of pixels

        

    #now sort list according to slice location
    imagesandSlices.sort(key= lambda x: x[1])
    for idx, item in enumerate(imagesandSlices):
        array[0,idx,:,:] = item[0]
        array[3,idx,:,:] = item[1]
        
        for x_idx in range(xLen):
            for y_idx in range(yLen):
                x = ipp[0] + x_idx*pixelSpacing[0]
                y = ipp[0] + y_idx*pixelSpacing[1]
                array[1,idx,x_idx,y_idx] = x
                array[2,idx,x_idx,y_idx] = y
        print("Collected PET SUV data for " + str(idx+1) + "/" + str(len(imagesandSlices)) + " DICOM files")         
    return array
    print("Got SUV Array for " + patientName)


def FindStructure(metadata, organ, invalidStructures = []):
    #Here we take the string for the desired structure (organ) and find the matching structure for each patient. 
    #The algorithm is to first make sure that the organ has a substring matching the ROI with at least 3 characters,
    #then out of the remaining possiblities, find top 3 closest fits with damerau levenshtein algorithm, then check to make sure that they are allowed to match according to rules defined in AllowedToMatch(). There should then ideally
    # be only one possible match, but if there are two, take the first in the list.   
    #Get a list of all structures in structure set: 
    unfilteredStructures = []
    for element in metadata:
        if element.get("ROIName").lower() not in invalidStructures:
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
    closestStrings = [["No Match",100],["No Match",100],["No Match",100]] #has to be in the top 3 closest strings to check next conditions
    for structure in structures:
        closeness = StringDistance(structure, organ)
        closestStrings.sort(key=itemgetter(1)) #Sort by the closeness value, and not the structure names
        for value in range(len(closestStrings)):
            if closeness < closestStrings[value][1]: #If closer than a value already in the top 3
                closestStrings[value] = [structure, closeness]
                break
    
    if closestStrings[0][1] == 100:
        return "No Match", 1111    
    #Now return the organ that is remaining and has closest string
    for element in metadata:
        if element.get("ROIName").lower() == closestStrings[0][0]:
            roiNumber = element.get("ROINumber")
    try:
        return closestStrings[0][0], roiNumber 
    except:
        return "No Match", 1111 #error code for unfound match.    

def AllowedToMatch(s1, s2):
    s1 = s1.lower()
    s2 = s2.lower()
    allowed = True
    keywords = []
    #You can't have only one organ with one of these keywords...
    keywords.append("prv")
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

    #Cant have left and no l in other, or rightt and no r
    if "left" in s1:
        if "l" not in s2:
            allowed = False      
    if "left" in s2:
        if "l" not in s1:
            allowed = False    
    #its tricky matching up left and right organs sometimes with all the conventions used... this makes sure that both are left or both are right
    if ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3]):
        if not (("lpar" in s2) or ("lsub" in s2) or ("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("_lt" in s2) or ("-lt " in s2) or ("lt " == s2[0:3])):   
            allowed = False  
    if (("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("-lt " in s2)or ("lt " == s2[0:3])):  
        if not (("lpar" in s1) or ("lsub" in s1) or ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3])):
            allowed = False        
    
    if ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)or ("right" in s1):
        if not (("rpar" in s2) or ("rsub" in s2) or ("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("_rt" in s2) or ("-rt" in s2) ):   
            allowed = False
    if (("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("-rt" in s2) ): 
        if not (("rpar" in s1) or ("rsub" in s1) or ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)):
            allowed = False
    return allowed


def StringDistance(s1, s2):
    return damerauLevenshtein(s1,s2,similarity=False)

    return d[lenstr1-1,lenstr2-1]
def LongestSubstring(s1,s2):
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

def GetContourMasks(rp_contours, CT_Array):
    numImagesSlices, xLen, yLen = CT_Array.shape[1:]
    contourMasks = np.zeros((2,numImagesSlices,xLen,yLen)) #2 channels, one for filled and one for unfilled

    rp_contours = CartesianToPixelCoordinates(rp_contours, CT_Array)
    for idx in range(numImagesSlices):#loop through all slices creating a mask for the contours
        for contour in rp_contours:
            if int(round(contour[0][2], 1)*100) == int(round(CT_Array[3,idx,0,0], 1)*100): #if contour is on the current slice
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
    return contourMasks                


def CartesianToPixelCoordinates(contours, array):
    #convert x and y values for a contour into the pixel indices where they are on the pet array
    xVals = array[1,0,:,0]
    yVals = array[2,0,0,:]
    for contour in contours: 
        for point in contour:
            point[0] = min(range(len(xVals)), key=lambda i: abs(xVals[i]-point[0]))
            point[1] = min(range(len(yVals)), key=lambda i: abs(yVals[i]-point[1]))
    return contours        

def ImageUpsizer(array, factor):
    """Supersizes an array by the factor given.   

    Args:
        array (2D numpy array): an array to be supersized
        factor (int): the amount to supersize the array by 

    Returns:
        newArray (2D numpy array): the supersized array

    """

    #Take an array and supersize it by the factor given
    xLen, yLen = array.shape
    newArray = np.zeros((factor * xLen, factor * yLen))
    #first get the original values in to the grid: 
    for i in range(xLen):
        for j in range(yLen):
            newArray[i * factor, j * factor] = array[i,j]
    #sample along first dim
    for j in range(yLen):
        for i in range(xLen - 1):
            insert = 1 
            while insert <= factor - 1:
                newArray[i * factor + insert, j * factor] = newArray[i * factor, j * factor] + (insert / factor) * (newArray[(i+1) * factor, j * factor]- newArray[i * factor, j * factor])
                insert += 1
    #sample along second dim
    for i in range(xLen * factor):
        for j in range(yLen - 1):
            insert = 1 
            while insert <= factor - 1:
                newArray[i, j * factor + insert] = newArray[i, j * factor] + (insert / factor) * (newArray[i, (j+1) * factor]- newArray[i, j * factor])
                insert += 1
    return newArray
if __name__ == '__main__':
    GetContours("/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1")
    #GetSUVArray("/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1")