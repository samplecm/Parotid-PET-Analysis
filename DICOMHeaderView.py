import os 
import pydicom 

path = "/media/caleb/WDBlue/PET_PSMA/PSMA_Analysis/SG_PETRT/3/CT_143516/1.2.826.0.1.3680043.10.740.1011040719149058956488316607536350590.dcm"

data = pydicom.dcmread(path)

print(data)
