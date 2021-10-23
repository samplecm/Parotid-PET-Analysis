#!/usr/bin/env bash
clear
echo "Clearing Terminal"
#This script gets radiomic csvs for parotids and submandibulars using dicomautomaton
cd SG_PETRT 

for patient in */; do
    patientPath="$PWD/$patient"
    echo $patientPath
    if [[ "$patient" == *"Contours"* ]]; then
        echo "skipping $patient"
        continue
    fi
    echo " "
    echo "Preparing patient: $patient"
    cd $patient
    #first loop through and get the path to the struct file
    for file in *; do
        if [[ "$file" == *"STRUCT"* ]]; then
            structPath="$PWD/$file"
            structFile=$file
        fi
    done
    #Now copy this structfile into the pet folder 
    for directory in */; do
        if [[ "$directory" == *"PET"* ]]; then
            pet_dir="$PWD/$directory"
            cd $directory
            #copy struct file to pet folder
            cp -v $structPath "$pet_dir$structFile"
            echo "copied $structFile to the PET folder."  
            echo "Left Par"
            #Now prep dicomautomaton
            rm "lp_features.csv"
            rm "rp_features.csv"        
            rm "ls_features.csv"
            rm "rs_features.csv"
            dicomautomaton_dispatcher * \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=lp_features.csv \
            -p NormalizedROILabelRegex='left parotid'    
            echo "Right Par"
            echo "$patientPath/lp_features.csv"
            cp -v "lp_features.csv" "$patientPath/Radiomics/lp_features.csv"
            rm "lp_features.csv"  
            cd ..
            cd $directory
            dicomautomaton_dispatcher * \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=rp_features.csv \
            -p NormalizedROILabelRegex='right parotid'    
            echo "Left Sub"
            cp -v "rp_features.csv" "$patientPath/Radiomics/rp_features.csv"
            rm "rp_features.csv"  
            cd ..
            cd $directory
            dicomautomaton_dispatcher * \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=ls_features.csv \
            -p NormalizedROILabelRegex='left submandibular'  
            echo "Right Sub"
            cp -v "ls_features.csv" "$patientPath/Radiomics/ls_features.csv"
            rm "ls_features.csv"  
            cd ..
            cd $directory
            dicomautomaton_dispatcher * \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=rs_features.csv \
            -p NormalizedROILabelRegex='right submandibular' 
            cp -v "rs_features.csv" "$patientPath/Radiomics/rs_features.csv"
            rm "rs_features.csv"    
               

            #delete struct file from pet folder
            rm "$pet_dir$structFile"
            cd .. 
        fi
    done

    cd .. 
    
done
