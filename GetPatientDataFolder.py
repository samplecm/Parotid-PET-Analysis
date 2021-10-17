#This script goes through and creates a folder containing a directory for each patient that only has pet data and structures.
import os 
import shutil

newPath = os.path.join(os.getcwd(), "resortedPatients")
if not os.path.isdir(newPath):
    os.mkdir(newPath)

patientDir = os.path.join(os.getcwd(), "SG_PETRT")
patientFiles = os.listdir(patientDir)
for patient in patientFiles:
    patientPath = os.path.join(newPath, patient)
    if not os.path.isdir(patientPath):
        os.mkdir(patientPath)

    dirs = os.listdir(os.path.join(patientDir, patient))
    #Get the pet folder and structure file
    for dir in dirs:
        if "PET" in dir:
            petPath = os.path.join(patientDir, patient, dir)
            for petImage in os.listdir(petPath):
                shutil.copyfile(os.path.join(petPath, petImage), os.path.join(patientPath, petImage))
        elif "RTSTRUCT" in dir:
            rtPath = os.path.join(patientDir, patient, dir)
            shutil.copyfile(rtPath, os.path.join(patientPath, dir))