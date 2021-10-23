import os 
import pydicom 

path = "/media/caleb/WDBlue/PET_PSMA/PSMA_Analysis/SG_PETRT/1/RTSTRUCT_SG-03"

data = pydicom.dcmread(path)

print(data)


# path = "/media/caleb/WDBlue/Organogenesis/Organogenesis/Patient_Files/P4/CT_000140.dcm"

# data = pydicom.dcmread(path)

# print(data)

