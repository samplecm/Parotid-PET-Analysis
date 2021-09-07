import GetImageData
from GetImageData import GetParotidSUVAnalysis, GetPatientName
import Segmentation
import Visuals
from Patient import Patient
from pydicom import dcmread
import pickle
import os
import numpy as np

def GetPatient(patientPath, patientNum):
   print("Loading patient data")
   patientName = GetPatientName(patientPath)
   print("  Getting left parotid contours")
   lp_contours = GetImageData.GetContours(patientPath, "Left Parotid")
   print("  Getting right parotid contours")
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
      print("Upsizing the PET arrays to match CT resolution")
      for idx in range(PET_Array.shape[1]):
         PET_Array_SS[0,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[0,idx,:,:],[512,512])
         PET_Array_SS[1,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[1,idx,:,:],[512,512])
         PET_Array_SS[2,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[2,idx,:,:],[512,512])
         PET_Array_SS[3,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[3,idx,:,:],[512, 512])
         print("    " + str(idx+1) + "/" + str(PET_Array.shape[1]) + " upsized.")     
      patient.PETArray = PET_Array_SS
      with open(os.path.join(patientPath, "PatientData.txt"), "wb") as fp:
         pickle.dump(patient, fp)
   patient.LeftParotid = lp_contours
   patient.RightParotid = rp_contours

   patient.LeftParotidMasks = GetImageData.GetContourMasks(lp_contours.wholeROI.copy(), patient.PETArray)
   patient.RightParotidMasks = GetImageData.GetContourMasks(rp_contours.wholeROI.copy(), patient.PETArray)
         
   return patient
for i in range(17,31):

   patientPath = "/media/calebsample/Data/PET PSMA/PSMA Analysis/SG_PETRT/" + str(i)
   patient = GetPatient(patientPath, 1)
   #Visuals.plotStructure(patient.RightParotid.segmentedContours18)
   #Visuals.PlotPETwithParotids(patient)
   GetParotidSUVAnalysis(patient)
   #Visuals.PlotCTwithParotids(patient)
   #Visuals.PlotPETwithParotids(patient, plotSubSegs=False)
   #Visuals.plotStructure(patient.LeftParotid.segmentedContours18)

