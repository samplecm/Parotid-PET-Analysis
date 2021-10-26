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
    #Now get the subsegments struct file
    for file in *; do
        if [[ "$file" == *"Subsegs"* ]]; then
            subSegsPath="$PWD/$file"
            subSegsFile=$file
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
            #copy subsegs file to pet folder
            cp -v $subSegsPath "$pet_dir$subSegsFile" 
            echo "copied $subSegsFile to the PET folder."   

            #Now prep dicomautomaton
            
            dicomautomaton_dispatcher *STRUCT* \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=lp_features.csv \
            -p NormalizedROILabelRegex="left parotid"    
            cp -v "lp_features.csv" "${patientPath}Radiomics/lp_features.csv"
            rm "lp_features.csv"  
            cd ..
            cd $directory
            dicomautomaton_dispatcher *STRUCT* \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=rp_features.csv \
            -p NormalizedROILabelRegex="right parotid"    
            echo "Left Sub"
            cp -v "rp_features.csv" "${patientPath}Radiomics/rp_features.csv"
            rm "rp_features.csv"  
            cd ..
            cd $directory
            dicomautomaton_dispatcher *STRUCT* \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=ls_features.csv \
            -p NormalizedROILabelRegex='left submandibular'  
            echo "Right Sub"
            cp -v "ls_features.csv" "${patientPath}Radiomics/ls_features.csv"
            rm "ls_features.csv"  
            cd ..
            cd $directory
            dicomautomaton_dispatcher *STRUCT* \
            -o ExtractRadiomicFeatures \
            -p FeaturesFileName=rs_features.csv \
            -p NormalizedROILabelRegex='right submandibular' 
            cp -v "rs_features.csv" "${patientPath}Radiomics/rs_features.csv"
            rm "rs_features.csv"    
            #delete struct file from pet folder
            rm "$pet_dir$structFile"   

            cd .. 

            #Now do this for all subsegments of both parotids
            for subseg_idx in {1..18}; do
                cd $directory
                csvNameRight="rParSub_$subseg_idx.csv"
                csvNameLeft="lParSub_$subseg_idx.csv"
                subNameRight="subseg_RPar$subseg_idx"
                subNameLeft="subseg_LPar$subseg_idx"

                dicomautomaton_dispatcher *Subseg* \
                -o ExtractRadiomicFeatures \
                -p FeaturesFileName=$csvNameRight \
                -p NormalizedROILabelRegex="$subNameRight" 
                cp -v "$csvNameRight" "$patientPath/Radiomics/$csvNameRight"
                rm "$csvNameRight"    
                cd ..

                cd $directory
                dicomautomaton_dispatcher *Subseg* \
                -o ExtractRadiomicFeatures \
                -p FeaturesFileName=$csvNameLeft \
                -p NormalizedROILabelRegex="$subNameLeft" 
                cp -v "$csvNameLeft" "$patientPath/Radiomics/$csvNameLeft"
                rm "$csvNameLeft"  
                cd ..



            done    
            cd $directory
            #delete subsegs file from pet folder
            rm "$pet_dir$subSegsFile"
            cd ..

        fi
    done

    cd .. 
    
done
