import GetImageData
from GetImageData import GetPatientName
from Data_Analyzing import ParotidSUVAnalysis, SubmandibularSUVAnalysis
import Segmentation
import Visuals
from Patient import Patient
from pydicom import dcmread
import pickle
import os
import numpy as np
import Data_Analyzing
from Data_Analyzing import DicomSaver

parentDirectory = os.getcwd()

def GetPatient(patientPath, patientNum, cleanProcess=False):
   print("Loading patient data")
   patientName = GetPatientName(patientPath)

   print("  Getting left parotid contours")
   lp_contours = GetImageData.GetContours(patientPath, "Left Parotid")

   print("  Getting right parotid contours")
   rp_contours = GetImageData.GetContours(patientPath, "Right Parotid")

   print("  Getting left submandibular contours")
   ls_contours = GetImageData.GetContours(patientPath, "Left Submandibular", subsegmentation=[1,1,1])

   print("  Getting right submandibular contours")
   rs_contours = GetImageData.GetContours(patientPath, "Right Submandibular", subsegmentation=[1,1,1])
   #Visuals.plotStructure(lp_contours.segmentedContours18)
   # Visuals.plotStructure(ls_contours.segmentedContours8)
   # Visuals.plotStructure(rs_contours.segmentedContours8)
   if cleanProcess == True:
      patient = Patient(patientName)
      CT_Array = GetImageData.GetCTArray(patientPath)
      PET_Array = GetImageData.GetSUVArray(patientPath)
      patient.path = patientPath 
      patient.patientNum = patientNum
      patient.CTArray = CT_Array
      PET_Array_SS = np.zeros((5, PET_Array.shape[1], 512, 512))
      print("Upsizing the PET arrays to match CT resolution")
      for idx in range(PET_Array.shape[1]):
         PET_Array_SS[0,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[0,idx,:,:],[512,512])
         PET_Array_SS[1,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[1,idx,:,:],[512,512])
         PET_Array_SS[2,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[2,idx,:,:],[512,512])
         PET_Array_SS[3,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[3,idx,:,:],[512,512])
         PET_Array_SS[4,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[4,idx,:,:],[512,512])
         print("    " + str(idx+1) + "/" + str(PET_Array.shape[1]) + " upsized.")     
      patient.PETArray = PET_Array_SS
   else:
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
         PET_Array_SS = np.zeros((5, PET_Array.shape[1], 512, 512))
         print("Upsizing the PET arrays to match CT resolution")
         for idx in range(PET_Array.shape[1]):
            PET_Array_SS[0,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[0,idx,:,:],[512,512])
            PET_Array_SS[1,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[1,idx,:,:],[512,512])
            PET_Array_SS[2,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[2,idx,:,:],[512,512])
            PET_Array_SS[3,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[3,idx,:,:],[512,512])
            PET_Array_SS[4,idx,:,:] = GetImageData.ImageUpsizer(PET_Array[4,idx,:,:],[512,512])
            print("    " + str(idx+1) + "/" + str(PET_Array.shape[1]) + " upsized.")     
         patient.PETArray = PET_Array_SS
      
   patient.LeftParotid = lp_contours
   patient.RightParotid = rp_contours
   patient.LeftSubmandibular = ls_contours
   patient.RightSubmandibular = rs_contours
   with open(os.path.join(patientPath, "PatientData.txt"), "wb") as fp:
         pickle.dump(patient, fp)
   DicomSaver(patientPath, ["parotids"], 18)

   patient.LeftParotidMasks = GetImageData.GetContourMasks(lp_contours.wholeROI.copy(), patient.PETArray)
   patient.RightParotidMasks = GetImageData.GetContourMasks(rp_contours.wholeROI.copy(), patient.PETArray)
   patient.LeftSubmandibularMasks = GetImageData.GetContourMasks(ls_contours.wholeROI.copy(), patient.PETArray)
   patient.RightSubmandibularMasks = GetImageData.GetContourMasks(rs_contours.wholeROI.copy(), patient.PETArray)
         
   return patient

# Data_Analyzing.Population_Metrics_SM()  
# Data_Analyzing.Population_Metrics_Parotid()  
Visuals.CorrelationPlot("submandibular")
# for i in range(10,31):
#    print("Loading Patient: " + str(i))
#    patientPath = os.path.join(parentDirectory, "SG_PETRT" , str(i))
#    patient = GetPatient(patientPath, 1, cleanProcess = False)

#    #Visuals.PlotSUVs(patient)
#    #Visuals.plotStructure(patient.RightParotid.segmentedContours18)
#    #Visuals.PlotPETwithParotids(patient)
#    SubmandibularSUVAnalysis(patient)
#    ParotidSUVAnalysis(patient)
   
   
#    #Visuals.PlotCTwithParotids(patient)
#    #Visuals.PlotPETwithParotids(patient, plotSubSegs=False)
#    #Visuals.plotStructure(patient.LeftParotid.segmentedContours18)
# Data_Analyzing.Population_Metrics_SM()   
# Data_Analyzing.Population_Metrics_Parotid()  

