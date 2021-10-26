#!/bin/bash
clear
echo "Clearing Terminal"
#This script gets radiomic csvs for parotids and submandibulars using dicomautomaton
cd SG_PETRT 

for patient in */; do
    echo $patient
    cd $patient
    #first loop through and get the path to the struct file
    python - << EOF
import os 
import pydicom
currentPath = os.getcwd()
for file in os.listdir(currentPath):
    if "STRUCT" in file:
        print(file)
        absFile = os.path.join(currentPath, file)
        data = pydicom.dcmread(absFile)
        name = data[0x0010,0x0010].value
        if name == "SG-30":
            print("Found patient " + str(name))   
EOF
cd ..
done