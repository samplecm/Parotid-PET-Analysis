
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import math
from GetImageData import CloneList
from numpy.lib.type_check import imag
import GetImageData
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from Patient import Patient
from PIL import Image
import copy

def plotStructure(structure):
    print("In PlotStructure")
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    colour_idx = 0
    colours = ['r', 'b', 'g', 'y', 'm', 'c', 'k']
    for substructure in structure:
        
        colour = colours[colour_idx]
        colour_idx = (colour_idx + 1) % 7
        for contour in substructure:
            
            x = []
            y = []
            z = []
            for point in contour:
                x.append(point[0])
                y.append(point[1])
                z.append(point[2])
            ax.plot(x,y,z, colour)    
        
    plt.show()
    print("")


def PlotSUVs(patient):
    #this function creates a slider plot of all suv slice images 
    global fig, ax, sliderT, sliderO, array, img

    #First get the suvs as a 3d numpy array. 
    array = GetImageData.GetSUVArray(patient)
    #this function takes a 4d array with the 3rd index being the slice and the 4th index being the b value 
    fig, ax = plt.subplots()
    img = ax.imshow(array[0,0,:,:], origin = 'upper')

    axT = fig.add_axes([0.2, 0.95, 0.65, 0.03])

    sliderT = Slider(axT, 'Slice', 0, array.shape[1]-1, valinit=0, valfmt='%i')


    sliderT.on_changed(update)

    plt.show()

def update(val):
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
    PlotSUVs("/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1")
                                    
