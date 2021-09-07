import pydicom 
import os
import pathlib
import glob
import CT
import PET
import matplotlib.pyplot as plt

patientNames = [] 
for patient in range(1,29):
    path = os.path.join(pathlib.Path(__file__).parent.absolute(), "SG_PETRT", str(patient))

    #Make a list to store CT and PET objects
    CTs = []
    PETs = []
    
    imageFolders = os.listdir(path)
    for imageFolder in imageFolders:
        files = sorted(glob.glob(os.path.join(path, imageFolder, "*")))
        for file in files: 
            if "PET" in file:
                data = pydicom.dcmread(file)
                patientName = data[0x0010,0x0020].value
                if patientName not in patientNames:
                    patientNames.append(patientName)

                #image = data.pixel_array
                #plt.imshow(data.pixel_array, cmap=plt.cm.bone)  
                #plt.show()
                #ct = CT(image, file, timeStamp, ipp, iop)
patientNames.sort()
print("finished")




