import math
import Contours

def OrganChopper(contours, numCuts, organName):
    print("In organChopper")
    #chop into 18ths
    numCutsX = numCuts[0]
    numCutsY = numCuts[1]
    numCutsZ = numCuts[2]

    ZChop(contours, numCutsZ) #Make axial cuts and store in contours.segmentedContours3

    contoursY = []
    if numCutsY > 0:
        for i in range(len(contours.segmentedContours3)):
            temp = YChop(contours.segmentedContours3[i], numCutsY)
            for j in range(len(temp)):
                contoursY.append(temp[j].copy())
    else:
        contoursY = contours.segmentedContours3.copy() 

    finalContours = []
    if numCutsX > 0:
        for i in range(len(contoursY)):
            temp = XChop(contoursY[i], numCutsX)
            for j in range(len(temp)):
                finalContours.append(temp[j].copy())
    else:
        finalContours = contoursY

    finalContours = ReOrder(finalContours, organName, numCuts)    

    if len(finalContours) == 18:
        contours.segmentedContours18 = finalContours

          


    #XChop(contours, numCutsX)


def ZChop(contours, numCutsZ):
    print("In ZChop")
    zCuts = BestCutZ(contours, numCutsZ) #Get the list of z values where structure is to be chopped

    newContoursList = []

    for i in range(len(zCuts)):
        contoursZ = ClosestContourZ(zCuts[i], contours) #Get the closest contour above and below the cut
        newContour = []
        for j in range(len(contours.wholeROI[contoursZ[0]])): #Go over all points in closest contour

            point1 = contours.wholeROI[contoursZ[0]][j]
            #Now get the closest point in second closest contour
            point2 = ClosestPoint(point1, contours.wholeROI[contoursZ[1]])
            #now interpolate between the two
            newPoint = InterpolateXY(point1, point2, zCuts[i])

            #add new point to new contour
            newContour.append(newPoint)
        #add this new contour to newContoursList
        newContoursList.append(newContour)    
    #Now add all the original contours to the newcontourslist and sort it by z value
    for i in range(len(contours.wholeROI)):
        newContoursList.append(contours.wholeROI[i])
    newContoursList.sort(key=GetZVal)

    #Now create a list of contour lists for each axial slice
    axialContours = []
    for i in range(numCutsZ + 1):
        axialContours.append([])

    for i in range(len(newContoursList)):
        for j in range(numCutsZ):
            if newContoursList[i][0][2] <= zCuts[j]: #z value of first point in contour
                if j == 0:
                    axialContours[j].append(newContoursList[i].copy())
                elif newContoursList[i][0][2] >= zCuts[j-1]: #Needs to be larger than the cut below to be in a middle or top region
                    axialContours[j].append(newContoursList[i].copy())

            if newContoursList[i][0][2] >= zCuts[-1]: #z value of first point in contour
                if j == 0:
                    axialContours[-1].append(newContoursList[i].copy())        

            

            
    contours.segmentedContours3 =  axialContours  

    return axialContours           
            

def GetZVal(e): 
    #function for returning z value of a contour
    if len(e) == 0:
        raise ValueError("Contour of 0 points has no z value")
    return e[0][2]

def GetContourArea(contour):
    area = 0
    if len(contour) == 0:
        return 0

    for i in range(len(contour)-1):
        area += contour[i][0] * contour[i+1][1] - contour[i+1][0] * contour[i][1]
    area += contour[-1][0] * contour[0][1] - contour[0][0] * contour[-1][1]    
    return 0.5 * abs(area)

def BestCutZ(contours, numCuts):
    print("In BestCutZ")
    zCuts = [] 
    numContours = len(contours.wholeROI)
    totalVolume = 0
    deltaZ = abs(contours.wholeROI[0][0][2]- contours.wholeROI[1][0][2]) #axial distance between slices 
    contourAreas = []
    for i in range(numContours):
        area = GetContourArea(contours.wholeROI[i])
        contourAreas.append(area)
        if i != numContours - 1:  #Dont include the final area in the volume calculation, since it would be over-extending
            totalVolume += area * deltaZ

    #Now find the correct cut spots which result in equal volume subsegments  
    for i in range(1, numCuts + 1):
        volumeGoal = i / (numCuts + 1)    #This is the ratio of volume that each subsegment should consist of (python integer division returns a float)
        contIndex = 0
        subVolume = 0
        while (subVolume < volumeGoal * totalVolume): #This loop will result in contIndex which tells us that the new slice should be somewhere between contourIndex and contourIndex - 1
            subVolume += contourAreas[contIndex] * deltaZ
            contIndex += 1
        volumeBelow = 0
        for j in range(contIndex - 1):
            avgArea = 0.5 * (contourAreas[j] + contourAreas[j+1])#first get the average area between two contours, used to approximate volume between the two (makes no difference while using equal slice diffs)
            volumeBelow += avgArea * deltaZ
        #now get the avg area for the slicing region: 
        avgArea = 0.5 * (contourAreas[contIndex - 1] + contourAreas[contIndex])

        zSlice = contours.wholeROI[contIndex - 1][0][2] + (volumeGoal * totalVolume - volumeBelow) / avgArea
        zCuts.append(zSlice)
    return zCuts    


def ClosestContourZ(z, contours):
    print("In ClosestContourZ")
    temp = 1000
    closestContours = [0,0]

    #First find closest contour
    for i in range(len(contours.wholeROI)):
        contourDistance = abs(contours.wholeROI[i][0][2] - z)    
        if contourDistance < temp:
            closestContours[0] = i
            temp = contourDistance


    #Now get second closest contour
    temp = 1000
    for i in range(len(contours.wholeROI)):
        if i == closestContours[0]: #Can't count the closest contour again
            continue        
        contourDistance = abs(contours.wholeROI[i][0][2] - z)

        if contourDistance < temp:
            closestContours[1] = i
            temp = contourDistance

    return closestContours        


def ClosestPoint(point, contour):
    m = 1000
    closestPointIdx = 1000

    for i in range(len(contour)):
        diff = math.sqrt((point[0] - contour[i][0])**2 + (point[1] - contour[i][1])**2)
        if diff < m:
            closestPointIdx = i
            m = diff
    if closestPointIdx == 1000:
        raise Exception("Could not find closest point.")
    return contour[closestPointIdx]    

def InterpolateXY(point1, point2, z):
    if point1[2] == z:
        return point1
    elif point2[2] == z: 
        return point2

    xSlope = (point2[0] - point1[0]) / (point2[2] - point1[2])
    ySlope = (point2[1] - point1[1]) / (point2[2] - point1[2])   

    newX = point1[0] + xSlope * (z - point1[2])
    newY = point1[1] + ySlope * (z - point1[2])

    newPoint = [newX, newY, z]

    return newPoint


def YChop(contours, numCutsY):  
    yCuts = BestCutY(contours, numCutsY)

    #add intersection points to contours
    for i in range(len(contours)):
        for j in range(len(yCuts)):
            contours[i] = AddIntersectionsY(contours[i], yCuts[j])
            #contours[i] = ClosedLooper(contours[i])

    #Now divide into separate parts
    finalContours = []
    divisions = [] #list for each y division for current contour

    #make the list the correct size so there is an item for each y div
    for div in range(len(yCuts)+ 1):
        finalContours.append([])

    for i in range(len(contours)):
        divisions.clear()
        #make the list the correct size
        for div in range(len(yCuts) + 1):
            divisions.append([])

        for y in range(len(yCuts)+1): #a section for each cut, +1
            for j in range(len(contours[i])):
                if y == 0:
                    if contours[i][j][1] <= yCuts[y]:
                        divisions[y].append(contours[i][j])
                elif y == len(yCuts):
                    if contours[i][j][1] >= yCuts[y-1]:
                        divisions[y].append(contours[i][j])
                else:
                    if contours[i][j][1] >= yCuts[y-1] and contours[i][j][1] <= yCuts[y]:
                        divisions[y].append(contours[i][j])
        for div in range(len(divisions)):
            #temp = ClosedLooper(divisions[div])
            finalContours[div].append(divisions[div].copy())
    return finalContours        

def BestCutY(contours, numCutsY):
    errorTolerance = 0.01
    area = 0
    maxY = -1000
    minY = 1000
    yCuts= []

    #get total area of contours
    for contour in contours:
        area += GetContourArea(contour)

    for cut in range(numCutsY):
        areaGoal = (cut + 1) / (numCutsY + 1)

        for j in range(len(contours)):
            if len(contours[j]) > 0:
                for row in range(len(contours[j])):
                    if contours[j][row][1] > maxY:
                        maxY = contours[j][row][1]
                    if contours[j][row][1] < minY:
                        minY = contours[j][row][1]

        tempContours = []
        cutContours = [] 

        error = 1000
        while error > errorTolerance:     
            tempContours.clear()
            yCut = (minY + maxY) / 2
            newArea = 0

            for i in range(len(contours)):
                cutContours.clear()
                temp = AddIntersectionsY(contours[i], yCut)
                #temp = ClosedLooper(temp)
                tempContours.append(temp)
                for j in range(len(tempContours[-1])):
                    if tempContours[-1][j][1] <= yCut:
                        cutContours.append(tempContours[i][j].copy())
                if len(cutContours) > 0:
                    #cutContours = ClosedLooper(cutContours)
                    newArea = newArea + GetContourArea(cutContours)
            error = abs((newArea / area) - areaGoal)   
            if (newArea / area < areaGoal):
                minY = yCut
            elif (newArea / area > areaGoal):
                maxY = yCut

             
        yCuts.append(yCut) 

        return yCuts

def ReOrder(contours, organName, numCuts):
    #Reorder from inferior --> superior, medial --> lateral, anterior --> posterior
    j = 0
    finalContours = []

    if 'l' in organName.lower():
        for i in range(len(contours)):
            if i % (numCuts[2] + 1) == 0 and i > 0:
                j = j + 1
            index = j % ((numCuts[0] + 1) * (numCuts[1] + 1) * (numCuts[2] + 1))     
            finalContours.append(contours[index])
            j = j + (numCuts[0] + 1) * (numCuts[1] + 1)
    else: 
        j = numCuts[0]
        for i in range(len(contours)):
            if j >= (numCuts[0] + 1) * (numCuts[1] + 1) * (numCuts[2] + 1):
                j-= (numCuts[0] + 1) * (numCuts[1] + 1) * (numCuts[2] + 1) 
                if j == -1:
                    j += (numCuts[0] + 1) * (numCuts[1] + 1)
            finalContours.append(contours[j])
            j = j + (numCuts[0]+1) * (numCuts[1] + 1)
    return finalContours

def AddIntersectionsY(contour, yCut):
    finalContour = contour.copy()
    numAdded =1 #start at one, increment after adding each point, to keep track of where to add additional point (add to index)
    z = contour[0][2]
    #index 0 outside of loop:
    if contour[0][1]  > yCut and contour[-1][1] < yCut:
        if contour[0][0] == contour[-1][0]: #if xs are same don't need to interpolate
            xNew = contour[0][0]
        else:        
            m = (contour[0][1]- contour[-1][1]) / (contour[0][0] - contour[-1][0])
            xNew = (yCut - contour[-1][1]) / m + contour[-1][0]
        finalContour = AddPoint(finalContour, 0, [xNew, yCut, z])
        numAdded = numAdded + 1
    if contour[0][1] < yCut and contour[-1][1] > yCut:
        if contour[0][0] == contour[-1][0]: #if xs are same don't need to interpolate
            xNew = contour[0][0]
        else:  
            m = (contour[-1][1] - contour[0][1]) / (contour[-1][0] - contour[0][0])
            xNew = (yCut - contour[0][1])/m + contour[0][0]
        finalContour = AddPoint(finalContour, len(contour), [xNew, yCut, z])    

    for i in range(0, len(contour)-1):
        if contour[i][1] < yCut and contour[i+1][1] > yCut:
            if contour[i][0] == contour[i+1][0]: #if xs are same don't need to interpolate
                xNew = contour[i][0]
            else:  
                m = (contour[i+1][1] - contour[i][1]) / (contour[i+1][0] - contour[i][0])
                xNew = (yCut - contour[i][1]) / m + contour[i][0]
            finalContour = AddPoint(finalContour, i + numAdded, [xNew, yCut, z])
            numAdded = numAdded + 1
        elif contour[i][1] > yCut and contour[i+1][1] < yCut:
            if contour[i][0] == contour[i+1][0]: #if xs are same don't need to interpolate
                xNew = contour[i][0]
            else:    
                m = (contour[i+1][1] - contour[i][1])/(contour[i+1][0] - contour[i][0])    
                xNew = (yCut - contour[i][1])/m + contour[i][0]
            finalContour = AddPoint(finalContour, i + numAdded, [xNew, yCut, z])
    return finalContour        

def AddPoint(contour, index, point):
    b = []
    if index > 0:
        for j in range(index):
            b.append(contour[j].copy())
    b.append(point)
    if index == len(contour):
        return b
    for j in range(index, len(contour)):
        b.append(contour[j].copy())
    return b 

def ClosedLooper(contour):
    if len(contour) == 0:
        return []
    numPoints = len(contour)
    x1 = contour[0][0]
    x2 = contour[numPoints-1][0]
    y1 = contour[0][1]
    y2 = contour[numPoints-1][1]

    if x1 != x2 or y1 != y2:
        contour.append(contour[0].copy())
    return contour    


def XChop(contours, numCuts):
    xCuts = BestCutX(contours, numCuts)
 
    for i in range(len(contours)): #First add intersection points
        for j in range(len(xCuts)):
            contours[i] = AddIntersectionsX(contours[i], xCuts[j])
        #contours[i] = ClosedLooper(contours[i])

        #now divide into separate parts
        finalContours = []
    divisions = [] #list for each y division for current contour

    #make the list the correct size so there is an item for each x div
    for div in range(len(xCuts)+ 1):
        finalContours.append([])

    for i in range(len(contours)):
        divisions.clear()
        #make the list the correct size
        for div in range(len(xCuts) + 1):
            divisions.append([])

        for x in range(len(xCuts)+1): #a section for each cut, +1
            for j in range(len(contours[i])):
                if x == 0:
                    if contours[i][j][0] <= xCuts[0]:
                        divisions[x].append(contours[i][j].copy())
                elif x == len(xCuts):
                    if contours[i][j][0] >= xCuts[x-1]:
                        divisions[x].append(contours[i][j].copy())
                else:
                    if contours[i][j][0] >= xCuts[x-1] and contours[i][j][0] <= xCuts[x]:
                        divisions[x].append(contours[i][j].copy())
        for div in range(len(divisions)):
            temp = ClosedLooper(divisions[div])
            finalContours[div].append(temp.copy())
    return finalContours   


def BestCutX(contours, numCuts):
    errorTolerance = 5e-5
    volume = 0
    xCuts= []

    #get total area of contours
    for i in range(len(contours)-1):
        if len(contours[i]) > 0 and len(contours[i+1]) > 0:
            deltaZ = abs(contours[i+1][0][2]-contours[i][0][2])
            volume = volume + (GetContourArea(contours[i]) + GetContourArea(contours[i+1])) * deltaZ

    for cut in range(numCuts):
        maxX = -1000
        minX = 1000
        areaGoal = (cut + 1) / (numCuts + 1)

        for j in range(len(contours)):
            if len(contours[j]) > 0:
                for row in range(len(contours[j])):
                    if contours[j][row][0] > maxX:
                        maxX = contours[j][row][0]
                    if contours[j][row][0] < minX:
                        minX = contours[j][row][0]

        tempContours = []
        cutContours = [] 
        numIters = 0
        error = 1000
        while error > errorTolerance and numIters < 1000:     
            
            tempContours.clear()
            xCut = (minX + maxX) / 2
            newAreas = []
            newVolume = 0

            for i in range(len(contours)):
                cutContours.clear()
                temp = AddIntersectionsX(contours[i], xCut)
                #temp = ClosedLooper(temp)
                tempContours.append(temp)
                for j in range(len(tempContours[-1])):
                    if tempContours[-1][j][0] <= xCut:
                        cutContours.append(tempContours[i][j].copy())
                if len(cutContours) > 0:
                    #cutContours = ClosedLooper(cutContours)
                    if len(cutContours) > 0:
                        newAreas.append([GetContourArea(cutContours), cutContours[0][2]]) #keep the z value as well for computing volume
                    else: 
                        newAreas.append([])    
            for area_idx in range(len(newAreas)-1):    
                if len(newAreas[area_idx]) > 0 and len(newAreas[area_idx+1]) > 0:
                    deltaZ = abs(newAreas[area_idx][1]-newAreas[area_idx+1][1])
                    newVolume = newVolume + (newAreas[area_idx][0] + newAreas[area_idx+1][0]) * deltaZ    
            error = abs((newVolume / volume) - areaGoal)   
            if (newVolume / volume < areaGoal): 
                minX = xCut            
                
            elif (newVolume / volume > areaGoal):
                maxX = xCut    
            #Re-adjust the max/min if stuck
               
            numIters = numIters + 1
  
        xCuts.append(xCut) 

    return xCuts

def AddIntersectionsX(contour, xCut):
    if len(contour) == 0:
        return []
    finalContour = contour.copy()
    numAdded =1 #start at one, increment after adding each point, to keep track of where to add additional point (add to index)
    z = contour[0][2]
    #index 0 outside of loop:
    if contour[0][0]  > xCut and contour[-1][0] < xCut:
        if contour[-1][1] == contour[0][1]: #if ys are same don't need to interpolate
            yNew = contour[0][1]
        else:        
            m = (contour[0][0]- contour[-1][0]) / (contour[0][1] - contour[-1][1])
            yNew = (xCut - contour[-1][0]) / m + contour[-1][1]
        finalContour = AddPoint(finalContour, 0, [xCut, yNew, z])
        numAdded = numAdded + 1
    if contour[0][0] < xCut and contour[-1][0] > xCut:
        if contour[0][1] == contour[-1][1]: #if ys are same don't need to interpolate
            yNew = contour[0][1]
        else:        
            m = (contour[-1][0] - contour[0][0]) / (contour[-1][1] - contour[0][1])
            yNew = (xCut - contour[0][0])/m + contour[0][1]
        finalContour = AddPoint(finalContour, len(contour), [xCut, yNew, z])    

    for i in range(0, len(contour)-1):
        if contour[i][0] < xCut and contour[i+1][0] > xCut:
            if contour[i+1][1] == contour[i][1]: #if ys are same don't need to interpolate
                yNew = contour[i][1]
            else:    
                m = (contour[i+1][0] - contour[i][0])/(contour[i+1][1] - contour[i][1])    
                yNew = (xCut - contour[i][0])/m + contour[i][1]
            finalContour = AddPoint(finalContour, i + numAdded, [xCut, yNew, z])
            numAdded = numAdded + 1
        elif contour[i][0] > xCut and contour[i+1][0] < xCut:
            if contour[i+1][1] == contour[i][1]: #if ys are same don't need to interpolate
                yNew = contour[i][1]
            else:    
                m = (contour[i+1][0] - contour[i][0])/(contour[i+1][1] - contour[i][1])    
                yNew = (xCut - contour[i][0])/m + contour[i][1]
            finalContour = AddPoint(finalContour, i + numAdded, [xCut, yNew, z])
            numAdded = numAdded + 1
    return finalContour        




                           








     




