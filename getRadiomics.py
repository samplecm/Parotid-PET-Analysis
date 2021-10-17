import os
import shutil

#Need to loop through patient files and perform radiomics using the structure files and pet images
patientFiles = os.path.join(os.getcwd(), "SG_PETRT")
patientFilesList = os.listdir(os.path.join(os.getcwd(), "SG_PETRT"))

for patient in patientFilesList: 
    dirs = os.listdir(os.path.join(patientFiles, patient))
    #Get the pet folder and structure file
    for dir in dirs:
        if "PET" in dir:
            petPath = os.path.join(patientFiles, patient, dir)
        elif "RTSTRUCT" in dir:
            rtPath = os.path.join(patientFiles, patient, dir)
            rtName = dir
    
    #put the structure file in its own folder to pass as input folder
    contourFolder = os.path.join(patientFiles, patient, "ContourFolder")
    if not os.path.isdir(contourFolder):
        os.mkdir(contourFolder)

    shutil.copyfile(rtPath, os.path.join(contourFolder, rtName))
    #Now can get the radiomics and save in a new folder.
    newFolder = os.path.join(patientFiles, patient, "Radiomics")
    if not os.path.isdir(newFolder):
        os.mkdir(newFolder)
    slicerCommand = "/home/caleb/Programs/Slicer/Slicer --no-main-window --python-script /home/caleb/Programs/SlicerRT/BatchProcessing/BatchStructureSetConversion.py --input-folder " + contourFolder + " --output-folder " \
         +  newFolder + " --ref-dicom-folder " + petPath
    os.system(slicerCommand)



