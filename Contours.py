class Contours(object):
    segmentedContours3 = []
    segmentedContours8 = []
    segmentedContours18 = []
    segmentedContours96 = []
    
    def __init__(self, name, dicomName, contours):
        self.wholeROI = contours
        self.roiName = name
        self.dicomName = dicomName
        
