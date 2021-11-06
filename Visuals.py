
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import math
from Data_Analyzing import CloneList
from numpy.lib.type_check import imag
import GetImageData
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from Patient import Patient
from PIL import Image
import copy
from Data_Analyzing import GetListRank
import NormalizeImportance


def ClusterPlot(array, labels):
    array = np.abs(array)
    numFeatures = array.shape[0]
    fig, ax = plt.subplots()
    img = ax.imshow(array)
    ax.set_xticks(range(0,array.shape[0]))
    ax.set_xticklabels(labels, rotation='vertical', fontsize=6)
    ax.set_yticks(range(0,array.shape[0]))
    ax.set_yticklabels(labels, fontsize=6)
    #ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
    plt.colorbar(img)
    plt.show()
    
    print("Created cluster plot.")


def plotStructure(structure):
    print("In PlotStructure")
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    minX = 1000
    minY = 1000
    minZ = 1000
    maxX = -1000
    maxY = -1000
    maxZ = -1000

    colour_idx = 0
    colours = ['r', 'b', 'g', 'y', 'm', 'c', 'k']
    for i in range(2,3):#len(structure)):
        substructure = structure[i]
        colour = colours[colour_idx]
        colour_idx = (colour_idx + 1) % 7
        for contour in substructure:          
            x = []
            y = []
            z = []
            for point in contour:
                if point[0] > maxX:
                    maxX = point[0]
                elif point[0] < minX:
                    minX = point[0]
                if point[1] > maxY:
                    maxY = point[1]
                if point[1] < minY:
                    minY = point[1]
                if point[2] > maxZ:
                    maxZ = point[2]
                if point[2] < minZ:
                    minZ = point[2]                     
                x.append(point[0])
                y.append(point[1])
                z.append(point[2])
                print(point)
            ax.plot(x,y,z, colour)    
    ax.set_xlim((minX-5, maxX+5))    
    ax.set_ylim((minY-5, maxY+5))    
    ax.set_zlim((minZ-5, maxZ+5))    

        
    plt.show()
    print("")


def PlotSUVs(patient : Patient):
    #this function creates a slider plot of all suv slice images 
    global fig, ax, sliderT, sliderO, array, img

    array = patient.PETArray
    
    #this function takes a 4d array with the 3rd index being the slice and the 4th index being the b value 
    fig, ax = plt.subplots()
    img = ax.imshow(array[0,0,:,:], origin = 'upper')

    axT = fig.add_axes([0.2, 0.95, 0.65, 0.03])

    sliderT = Slider(axT, 'Slice', 0, array.shape[1]-1, valinit=0, valfmt='%i')
    sliderV = Slider(axT, 'Slice', 0, array.shape[1]-1, valinit=0, valfmt='%i')


    sliderT.on_changed(updateSUVPlot)

    plt.show()

def updateSUVPlot(val):
    i = int(sliderT.val)
    im = array[0,i,:,:]

    img.set_data(im)
    # min = np.nanmin(array[:,:,i,j])
    # max = np.nanmax(array[:,:,i,j])
    # img.set_clim([min,max])


    #Only want to reload the colour bar for special images like D, D*, f. So we must clear the axis and make a new image for these types of images. Then when we go back to regular signal images need to recreate
    #specialImageRange = range(num_bVals - 3, num_bVals - 1)

    fig.canvas.draw()

    
def PlotPETwithParotids(patient : Patient, plotSubSegs=False):
    global fig, ax, sliderT, sliderO, combinedArray, rPar, lPar, img
    petArray = patient.PETArray
    #normalize the CT array to have minimum 0 and maximum 0.8:
    normalizedPetArray = NormalizeArray(CloneList(petArray), 0.8)
    
    if plotSubSegs:
        rPar = CloneList(patient.RightParotid.segmentedContours18)
        lPar = CloneList(patient.LeftParotid.segmentedContours18)
        lParMasks = GetImageData.GetContourMasks(CloneList(lPar[0]), petArray)
        rParMasks = GetImageData.GetContourMasks(CloneList(rPar[0]), petArray)
        combinedArray = MaskOnImage(normalizedPetArray.copy(), [lParMasks, rParMasks])
        for i in range(1,18):
            lParMasks = GetImageData.GetContourMasks(CloneList(lPar[i]), petArray)
            rParMasks = GetImageData.GetContourMasks(CloneList(rPar[i]), petArray)
            combinedArray = MaskOnImage(combinedArray, [lParMasks, rParMasks])
    else:        
        rPar = patient.RightParotidMasks
        lPar = patient.LeftParotidMasks
        combinedArray = MaskOnImage(normalizedPetArray.copy(),[lPar, rPar])
    #this function creates a slider plot of all suv slice images 
    
    #First get the suvs as a 3d numpy array. 

    #this function takes a 4d array with the 3rd index being the slice and the 4th index being the b value 
    fig, ax = plt.subplots()
    img = ax.imshow(combinedArray[0,0,:,:], origin = 'upper', vmin = 0, vmax = 1)

    axT = fig.add_axes([0.2, 0.95, 0.65, 0.03])
    sliderT = Slider(axT, 'Slice', 0, combinedArray.shape[1]-1, valinit=0, valfmt='%i')
    sliderT.on_changed(update)
    plt.show()

def update(val):
    i = int(sliderT.val)
    im = combinedArray[0,i,:,:]
    img.set_data(im)
    # min = np.nanmin(array[:,:,i,j])
    # max = np.nanmax(array[:,:,i,j])
    # img.set_clim([min,max])

    fig.canvas.draw()


    #Only want to reload the colour bar for special images like D, D*, f. So we must clear the axis and make a new image for these types of images. Then when we go back to regular signal images need to recreate
    #specialImageRange = range(num_bVals - 3, num_bVals - 1)

    fig.canvas.draw()

def PlotCTwithParotids(patient):
    global fig, ax, sliderT, sliderO, combinedArray, rPar, lPar, img
    ctArray = patient.CTArray
    #normalize the CT array to have minimum 0 and maximum 0.8:
    ctArray = NormalizeArray(ctArray, 0.8)
    print([np.amin(ctArray), np.amax(ctArray)])
    rPar = patient.RightParotidMasks
    lPar = patient.LeftParotidMasks
    combinedArray = MaskOnImage(ctArray.copy(),[lPar, rPar])
    #this function creates a slider plot of all suv slice images 
    
    #First get the suvs as a 3d numpy array. 

    #this function takes a 4d array with the 3rd index being the slice and the 4th index being the b value 
    fig, ax = plt.subplots()
    img = ax.imshow(combinedArray[0,0,:,:], origin = 'upper', vmin = 0, vmax = 1)

    axT = fig.add_axes([0.2, 0.95, 0.65, 0.03])
    sliderT = Slider(axT, 'Slice', 0, combinedArray.shape[1]-1, valinit=0, valfmt='%i')
    sliderT.on_changed(update)
    plt.show()

def update(val):
    i = int(sliderT.val)
    im = combinedArray[0,i,:,:]
    img.set_data(im)
    # min = np.nanmin(array[:,:,i,j])
    # max = np.nanmax(array[:,:,i,j])
    # img.set_clim([min,max])

    fig.canvas.draw()

def NormalizeArray(array, max):
    array_max = np.amax(array)
    array_min = np.amin(array)
    for i in range(array.shape[1]):
        array[:,i,:,:] = max * (array[:,i,:,:] - array_min) / (array_max - array_min)
    return array

def MaskOnImage(array, masks):
    newArray = array.copy()
    for i in range(array.shape[1]):
        image = newArray[0,i,:,:]
        for mask in masks:
            maskImage = mask[1,i,:,:]
            image = image + maskImage

        newArray[0,i,:,:] = image

    return newArray       

def CorrelationPlot(organ):
    if organ == "parotid":
        importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
            1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
            0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
            0.272984146341463,   0.0484897103658537,  0.035493902439024]
        leftSUVs = [10809, 13002, 12140, 13943, 16767, 16564, 16469, 18995, 17250, 13983, 16645, 15153, 17303, 20449, 18917, 18020, 21214, 19283]
        rightSUVs = [11820, 13593, 12279, 14884, 17103, 16998, 16550, 18245, 16755, 13209, 15898, 15443, 16076, 19399, 18955, 16905, 20393, 19244] 
    elif organ == "submandibular":
        importanceVals = NormalizeImportance.SubImportance()
        leftSUVs = [10.76,10.63,12.38,11.28,11.28,12.25,10.68,13.07,11.36]
        rightSUVs = [11.07,10.91,11.68,11.28,12.84,10.68,12.68,11.1]                
    leftSUVs = list(zip(leftSUVs, importanceVals))
    rightSUVs = list(zip(rightSUVs, importanceVals))

    #sort the lists by importance: 
    leftSUVs.sort(key=lambda x: x[1])
    rightSUVs.sort(key=lambda x: x[1])

    suvs_left = []
    suvs_right = []
    importanceVals.sort()
    for idx in range(len(leftSUVs)):
        suvs_left.append(leftSUVs[idx][0])
        suvs_right.append(rightSUVs[idx][0])
      
    fig = plt.figure()
    axLeft = fig.add_subplot(211)
    axLeft.scatter(importanceVals, suvs_left)
    axLeft.set_xlabel('Subsegment Importance')
    axLeft.set_ylabel('Subsegment Average SUV')
    if organ == "parotid":
        axLeft.set_title("Left Parotid: \n Spearman's Rank: -0.56 (p < 0.008) \n Pearson's Rank: -0.52 (p < 0.014)")
    else:
        axLeft.set_title("Left Sub: \n Spearman's Rank: -0.31 (p < 0.22) \n Pearson's Rank: -0.29 (p < 0.24)")    
    axRight = fig.add_subplot(212)
    axRight.scatter(importanceVals, suvs_right)
    axRight.set_xlabel('Subsegment Importance')
    axRight.set_ylabel('Subsegment Average SUV')
    if organ == "parotid":
        axRight.set_title("Right Parotid: \n Spearman's Rank: -0.512 (p < 0.015) \n Pearson's Rank: -0.477 (p < 0.023)")
    else:  
        axRight.set_title("Right Sub: \n Spearman's Rank: -0.24 (p < 0.28) \n Pearson's Rank: -0.29 (p < 0.24)")   
    plt.subplots_adjust(hspace=0.9)  
    plt.show()    




    



if __name__ == "__main__":
    # contour = []
    # radius = 1
    # offset = 0
    # for sub in range(3):
    #     contour.append([])
    #     for z in range(0,10):
    #         layer = []
    #         for i in range(1,22):
    #             angle = i/20 * 2 * math.pi 
    #             x = math.cos(angle) + offset
    #             y = math.sin(angle) + offset
    #             layer.append([x,y,z])
    #         contour[sub].append(layer)    
    #         radius = radius + 1
    #     offset = offset + 2.5

    # plotStructure(contour)    
    #PlotSUVs("/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1")
    CorrelationPlot()
                                    
