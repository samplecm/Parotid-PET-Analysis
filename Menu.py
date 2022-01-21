import GetImageData
from GetImageData import GetPatientName
from Data_Analyzing import KFold_Validation, ParotidSUVAnalysis, SubmandibularSUVAnalysis
import Segmentation
import Visuals
from Patient import Patient
from pydicom import dcmread
import pickle
import os
import numpy as np
import Data_Analyzing
from Data_Analyzing import DicomSaver
import Radiomics


parentDirectory = os.getcwd()

def GetPatient(patientPath, patientNum, cleanProcess=False, save=True):
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
         print("Trying to load previously processed patient data")
         with open(os.path.join(patientPath, "PatientData.txt"), "rb") as fp:
            patient = pickle.load(fp)
         PET_Array = patient.PETArray
         CT_Array = patient.CTArray
         print("Patient data loaded successfully")
         
      except:      
         print("Unable to load saved patient data. Processing DICOM files.")
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
   Visuals.plotSubsegments(patient.LeftParotid.segmentedContours18)
   with open(os.path.join(patientPath, "PatientData.txt"), "wb") as fp:
         pickle.dump(patient, fp)
   if save==True:      
      DicomSaver(patientPath, ["parotids", "submandibular"])

   patient.LeftParotidMasks = GetImageData.GetContourMasks(lp_contours.wholeROI.copy(), patient.PETArray)
   patient.RightParotidMasks = GetImageData.GetContourMasks(rp_contours.wholeROI.copy(), patient.PETArray)
   patient.LeftSubmandibularMasks = GetImageData.GetContourMasks(ls_contours.wholeROI.copy(), patient.PETArray)
   patient.RightSubmandibularMasks = GetImageData.GetContourMasks(rs_contours.wholeROI.copy(), patient.PETArray)
         
   return patient



#
   #Visuals.MeshPlot(patient.LeftParotid)
   # Visuals.PlotSUVs(patient)
   #Visuals.plotSubsegments(patient.LeftParotid.segmentedContours18)
   # Visuals.PlotPETwithParotids(patient)
   # SubmandibularSUVAnalysis(patient)
   # ParotidSUVAnalysis(patient, stat="95")
   # ParotidSUVAnalysis(patient, stat="5")
   
   
   #Visuals.PlotCTwithParotids(patient)
   #Visuals.PlotPETwithParotids(patient, plotSubSegs=True)
   #Visuals.plotStructure(patient.LeftParotid.segmentedContours18)
   
# labels_filtered = Radiomics.FeatureSelection("parotid")
# population_radiomics_stats = Radiomics.Get_Subseg_Features(labels_filtered, "parotid")
# Data_Analyzing.Model_Selection_aic(population_radiomics_stats, labels_filtered)

# importanceModel = Data_Analyzing.Model(population_radiomics_stats, degree=2)
# KFold_Validation(population_radiomics_stats, degree=10)
# Visuals.ModelScore_vs_Degree_Plot(population_radiomics_stats)
# prediction = importanceModel.predict(population_radiomics_stats[1,3,:])

# GetImageData.GetAvgAge()
# GetImageData.GetAvgWeight()
# GetImageData.GetSexes()

#Data_Analyzing.Population_Metrics_SM()   
#Data_Analyzing.Population_Metrics_Parotid(stat="mean")  
#Data_Analyzing.Population_Metrics_Parotid(stat="95")
#Data_Analyzing.Population_Metrics_Parotid(stat="5")
Visuals.CorrelationPlot("parotid")

print("finished")
