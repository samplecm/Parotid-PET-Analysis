import GetImageData
from GetImageData import GetParotidSUVs, GetPatientName
import Segmentation
import Visuals
from Patient import Patient
from pydicom import dcmread
import pickle
import os
import numpy as np

def GetPatient(patientPath, patientNum):
   patientName = GetPatientName(patientPath)
   lp_contours = GetImageData.GetContours(patientPath, "Left Parotid")
   rp_contours = GetImageData.GetContours(patientPath, "Right Parotid")
   try:
      with open(os.path.join(patientPath, "PatientData.txt"), "rb") as fp:
         patient = pickle.load(fp)
         PET_Array = patient.PETArray
         CT_Array = patient.CTArray
   except:      
      patient = Patient(patientName)
      CT_Array = GetImageData.GetCTArray(patientPath)
      PET_Array = GetImageData.GetSUVArray(patientPath)
      patient.path = patientPath 
      patient.patientNum = patientNum
      patient.CTArray = CT_Array
      patient.PETArray = PET_Array
      PET_Array_SS = np.zeros((4, PET_Array.shape[1], 512, 512))
      for idx in range(PET_Array.shape[1]):
         PET_Array_SS[0,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[0,idx,:,:],2)
         PET_Array_SS[1,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[1,idx,:,:],2)
         PET_Array_SS[2,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[2,idx,:,:],2)
         PET_Array_SS[3,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[3,idx,:,:],2)
      patient.PETArray = PET_Array_SS
      with open(os.path.join(patientPath, "PatientData.txt"), "wb") as fp:
         pickle.dump(patient, fp)
   patient.LeftParotid = lp_contours
   patient.RightParotid = rp_contours

   patient.LeftParotidMasks = GetImageData.GetContourMasks(lp_contours.wholeROI.copy(), patient.PETArray)
   patient.RightParotidMasks = GetImageData.GetContourMasks(rp_contours.wholeROI.copy(), patient.PETArray)
         
   return patient

patientPath = "/home/calebsample/Documents/UBC/PET PSMA/PSMA Analysis/SG_PETRT/1"
patient = GetPatient(patientPath, 1)
Visuals.plotStructure(patient.RightParotid.segmentedContours18)
#Visuals.PlotPETwithParotids(patient)
GetParotidSUVs(patient)
#Visuals.PlotCTwithParotids(patient)
Visuals.PlotPETwithParotids(patient)
Visuals.plotStructure(patient.RightParotid.segmentedContours18)

