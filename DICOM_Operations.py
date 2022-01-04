import os 
import pydicom 
def ShowMeta(path):
    data = pydicom.dcmread(path)
    print(data)

def Add_FOR_IDS():
    #adds frame of reference IDs from original dicoms to the new ones. 
    planning_files_path = "/media/caleb/WDBlue/PET_PSMA/Planning_Files_2"
    patients_path = "/media/caleb/WDBlue/PET_PSMA/PSMA_Analysis/SG_PETRT"
    for patient in sorted(os.listdir(patients_path), reverse=True):
        patient_path = os.path.join(patients_path, patient)
        planning_file_path = os.path.join(planning_files_path, patient, "STRUCT.dcm")
        #now load both the structures (limbus original , and one with added gtvs)
        limbus_file = None
        gtv_file = None
        for file in os.listdir(patient_path):
            if "RTSTRUCT" in file:
                limbus_file = pydicom.dcmread(os.path.join(patient_path, file))
            if "planning" in file.lower():    
                gtv_file = pydicom.dcmread(os.path.join(patient_path, file))
        if limbus_file == None:
            raise Exception("Limbus file not found")
        if gtv_file == None:
            raise Exception("GTV file not found")


        #Now get the sequence from limbus that needs to be added to gtv
            
        for_sequence = limbus_file[0x3006, 0x0010]

        #Now add to gtv file
        gtv_file.add_new([0x3006, 0x0010], "SQ", for_sequence)
        
        #save it
        if not os.path.exists(os.path.join(planning_files_path, patient)):
            os.mkdir(os.path.join(planning_files_path, patient))
        gtv_file.save_as(planning_file_path)
        print("")



if __name__ == "__main__":
    path = "/media/caleb/WDBlue/PET_PSMA/PSMA_Analysis/SG_PETRT/1/Subsegs.dcm"
    #ShowMeta()   
    Add_FOR_IDS() 






# path = "/media/caleb/WDBlue/Organogenesis/Organogenesis/Patient_Files/P4/CT_000140.dcm"

# data = pydicom.dcmread(path)

# print(data)

