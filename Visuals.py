
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import math
from Contours import Contours
from Data_Analyzing import CloneList
from numpy.lib.type_check import imag
import GetImageData
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from Patient import Patient
from PIL import Image
import copy
from Data_Analyzing import GetListRank, KFold_Validation
import NormalizeImportance
import plotly.graph_objects as go
import Contour_Operations

def ModelScore_vs_Degree_Plot(radiomics_data):
    x = np.linspace(1, 5, 5, dtype=int)
    y = []
    for i in x:
        score = KFold_Validation(radiomics_data, k=9, degree=i)
        y.append(score)

    plt.figure()
    plt.scatter(x,y)
    plt.xlabel("Model Degree")
    plt.ylabel("Mean Absolute Error")
    plt.xticks([0,1,2,3,4,5])
    plt.show()
    print("Created model score vs degree plot")

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

def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])

def plotSubsegments(structure):
    print("In PlotStructure")
    fig = plt.figure()
    ax : plt.Axes = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)


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
    importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
            1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
            0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
            0.272984146341463,   0.0484897103658537,  0.035493902439024]
    suv_values = [10.62, 8.7, 7.68, 8.83, 11.2, 10.4, 10.5, 12.7, 11.2, 8.9, 11, 9.8, 11.06, 13.64, 12.4, 11.7, 14.2, 12.8] 
    suv_values = [i- min(suv_values) for i in suv_values]  
    suv_values = [i / max(suv_values) for i in suv_values] 

    for i in range(len(structure)):
        substructure = structure[i]
        # colour = colours[colour_idx]
        # colour_idx = (colour_idx + 1) % 7
        importance = importanceVals[i]
        colour = Get_Colormap_RGB(importance)
        for c, contour in enumerate(substructure):   
              
            if len(contour) == 0:
                continue    
            # if c % 2 != 0:
            #     continue
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

            #ax.plot(x,y,z, colour)  
            points = [list(zip(x,y,z))]
            poly = a3.art3d.Poly3DCollection(points)  
            poly.set_color(colour)
            poly.set_edgecolor('k')
            ax.add_collection3d(poly)
    #now add with suv mapping
    for i in range(len(structure)):
        substructure = structure[i]
        # colour = colours[colour_idx]
        # colour_idx = (colour_idx + 1) % 7
        importance = suv_values[i]
        colour = Get_Colormap_RGB(importance)
        for c, contour in enumerate(substructure):   
                
            if len(contour) == 0:
                continue    
            # if c % 2 != 0:
            #     continue
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
                x.append(point[0]+50)
                y.append(point[1])
                z.append(point[2])

            #ax.plot(x,y,z, colour)  
            points = [list(zip(x,y,z))]
            poly = a3.art3d.Poly3DCollection(points)  
            poly.set_color(colour)
            poly.set_edgecolor('k')
            ax.add_collection3d(poly)

    ax.set_xlim((minX-5, maxX+5))    
    ax.set_ylim((minY-5, maxY+5))    
    ax.set_zlim((minZ-5, maxZ+5))    

    set_axes_equal(ax)
    ax.grid(False) 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([]) 
    plt.axis('off')
    plt.show()
    print("")

def Get_Colormap_RGB(val):
    
    # val_blue = min(0.5,val)
    # blue = (1-2*val_blue)

    # val_green = abs(val-0.5)
    # green = (1-2*val_green)

    # val_red = 1 - max(0.5,val)
    # red = (1-2*val_red)
    # i reveresed the order by accident... so reverse val
    val = 1 - val

    if val < 0.25:
        red = 1
        blue = 0.1
        green = 0.9*(val/0.25) + 0.1
    elif val < 0.5:
        red = -0.9*((val-0.25)/0.25) + 1
        green = 1
        blue = 0.1
    elif val < 0.75:
        red = 0.1
        green = 1    
        blue = 0.9*((val-0.5)/0.25) + 0.1
    else:
        red = 0.1
        blue = 1
        green = -0.9*((val-0.75)/0.25) + 1


    return (red, green, blue, 1)
    
def plotStructure_unfilled(structure):
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
    for i in range(len(structure)):
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

def CorrelationPlot(organ, stat="mean"):
    if organ == "parotid":
        importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
            1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
            0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
            0.272984146341463,   0.0484897103658537,  0.035493902439024]
        leftSUVs = [6.83, 8.7, 7.68, 8.83, 11.2, 10.4, 10.5, 12.7, 11.2, 8.9, 11, 9.8, 11.06, 13.64, 12.4, 11.7, 14.2, 12.8] 
        rightSUVs = [7.49, 8.99, 7.83, 9.45, 11.3, 10.75, 10.64, 12.17, 10.84, 8.34, 10.44, 10.05, 10.18, 12.86, 12.34, 10.96, 13.55, 12.7] 
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
    importanceVals = np.array(importanceVals)
    suvs_left = np.array(suvs_left)
    suvs_right = np.array(suvs_right)  
    #get line of best fits.
    m_l, b_l = np.polyfit(importanceVals, suvs_left, 1)  
    m_r, b_r = np.polyfit(importanceVals, suvs_right, 1)  


    fig = plt.figure()
    axLeft = fig.add_subplot(211)
    axLeft.scatter(importanceVals, suvs_left, color='darkorange')
    axLeft.plot(importanceVals, importanceVals*m_l + b_l, color='k')
    axLeft.set_xlabel('Relative Importance')
    axLeft.set_ylabel('Mean SUVbw')
    if organ == "parotid":
        axLeft.set_title("Left Parotid: \n Spearman's Rank: -0.56 (p < 0.007)")
    else:
        axLeft.set_title("Left Sub: \n Spearman's Rank: -0.31 (p < 0.22) \n Pearson's Rank: -0.29 (p < 0.24)")    
    axRight = fig.add_subplot(212)
    axRight.scatter(importanceVals, suvs_right, color='r')
    axRight.plot(importanceVals, importanceVals*m_r + b_r, color='k')
    axRight.set_xlabel('Relative Importance')
    axRight.set_ylabel('Mean SUVbw')
    if organ == "parotid":
        axRight.set_title("Right Parotid: \n Spearman's Rank: -0.54 (p < 0.01)")
    else:  
        axRight.set_title("Right Sub: \n Spearman's Rank: -0.24 (p < 0.28) \n Pearson's Rank: -0.29 (p < 0.24)")   
    plt.subplots_adjust(hspace=0.9)  
    plt.show()    

def MeshPlot(contour : Contours):
    #This creates a plot showing parotid subsegments with colour scheme according to importance
    importanceVals = [0.751310670731707,  0.526618902439024,   0.386310975609756,
                1,   0.937500000000000,   0.169969512195122,   0.538871951219512 ,  0.318064024390244,   0.167751524390244,
                0.348320884146341,   0.00611608231707317, 0.0636128048780488,  0.764222560975610,   0.0481192835365854,  0.166463414634146,
                0.272984146341463,   0.0484897103658537,  0.035493902439024]
    subsegs = contour.segmentedContours18
    parotid = contour.wholeROI
    parotid = Contour_Operations.AddInterpolatedPoints(copy.deepcopy(parotid))
    x_parotid = []
    y_parotid = []
    z_parotid = []
    for s, slice in enumerate(parotid):
        for p, point in enumerate(slice):
            x_parotid.append(point[0])
            y_parotid.append(point[1])
            z_parotid.append(point[2])
    fig = go.Figure(data=go.Mesh3d(
        x=x_parotid, y=y_parotid, z=z_parotid, color='green', alphahull=0.1, opacity = 1))
    fig.show()
    print()            




    



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
                                    
