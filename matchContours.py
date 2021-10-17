import os 
import pydicom
import shutil


#This file will find dicom rt files corresponding to all patients in the patientsdirectory and copy the file to the patient folders. 
#dicom structure files are placed in the patients directory in a folder called contours 

patientsDirectory = os.path.join(os.getcwd(), "SG_PETRT")
contoursFolder = os.path.join(patientsDirectory, "Contours")

patients = os.listdir(patientsDirectory)

for patient in patients:
    if patient == "Contours":
        continue
    patientPath = os.path.join(patientsDirectory, patient)
    #now go into the PET folder to get the name
    patientImageFolders = os.listdir(patientPath)
    for ImageTypefolder in patientImageFolders:
        if "PET" in ImageTypefolder:
            PETImages = os.listdir(os.path.join(patientPath, ImageTypefolder))
            petImagePath = os.path.join(patientPath, ImageTypefolder, PETImages[0])#Take the first pet image to get the name
            data = pydicom.dcmread(petImagePath)
            patientName = data[0x0010,0x0020].value
    
    #Now search in the contours folder for a dicom file with the same patient name.
    names = []
    for file in os.listdir(contoursFolder):
        structureFile = os.path.join(contoursFolder, file)
        data = pydicom.dcmread(structureFile)
        name = data[0x0010,0x0020].value
        names.append(name)
        if name == patientName:
            os.rename(structureFile, os.path.join(contoursFolder, str("RTSTRUCT_" + name)))
            shutil.copy(os.path.join(contoursFolder, str("RTSTRUCT_" + name)), patientPath)
names = sorted(names)

    
    
print("finished")

