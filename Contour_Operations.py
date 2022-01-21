import copy 


def AddInterpolatedPoints(orig_contours):
    """This makes sure that each slice of contours has at least 100 points
    Args: 
        contours (list): the contour list for a single patient
    """
    contours = copy.deepcopy(orig_contours)
    #Now add in sufficient number of points 
    for contour_idx in range(len(contours)): 
        contour = contours[contour_idx]
        numPointsOrig = len(contour)
        numPoints = len(contour)
        if numPoints > 100:
            continue
        if numPoints < 4:
            contours[contour_idx] = []
            continue


        pointIncreaseFactor = 1
        while numPoints < 100:  
            numPoints = numPoints + numPointsOrig
            pointIncreaseFactor = pointIncreaseFactor + 1      
        increasedContour = []
        for point_idx in range(len(contour)-1):    
            increasedContour.append(contour[point_idx].copy())
            for extraPointNum in range(pointIncreaseFactor):
                scaleFactor = extraPointNum / (pointIncreaseFactor + 1)
                newX = (contour[point_idx+1][0] - contour[point_idx][0]) * scaleFactor + contour[point_idx][0]
                newY = (contour[point_idx+1][1] - contour[point_idx][1]) * scaleFactor + contour[point_idx][1]
                z = contour[point_idx][2]
                newPoint = [newX, newY, z]
                increasedContour.append(newPoint)
        #Now do it for the last point connected to the first
        for extraPointNum in range(pointIncreaseFactor):
            scaleFactor = extraPointNum / (pointIncreaseFactor + 1)
            newX = (contour[0][0] - contour[-1][0]) * scaleFactor + contour[-1][0]
            newY = (contour[0][1] - contour[-1][1]) * scaleFactor + contour[-1][1]
            z = contour[-1][2]
            newPoint = [newX, newY, z]
            increasedContour.append(newPoint)        
        contours[contour_idx] = increasedContour

    return contours 