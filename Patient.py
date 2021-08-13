from math import nan


class Patient():
    def __init__(self, patientName):
        self.name = patientName
        self.PETArray = nan
        self.CTArray = nan
        self.LeftParotid = nan
        self.RightParotid = nan
        self.LeftParotidMasks = nan
        self.RightParotidMasks = nan