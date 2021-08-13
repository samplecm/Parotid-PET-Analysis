#This class will old the necessary PET data needed. 

class PET():
    def __init__(self, pixelData, fileName, timeStamp, ipp, iop):
        self.pixelData = pixelData
        self.fileName = fileName
        self.timeStamp = timeStamp
        self.ipp = ipp
        self.iop = iop
        


