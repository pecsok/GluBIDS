########################################
**pyGluCEST** : Current Processing Steps
########################################
Below are descriptions of the scripts currently used in pyGluCEST, divided into sections that could potentially provide the basis for processing modules for GluBIDS.


*****
Phase 1: Setting up the data.
*****
#### 1. copy_rawdata.py

Hard-coded Inputs:
- paths; expects rawdata folder to be in <subid>/<sesid> format
- required_files variable (line 27) designates which dicoms to copy over to new folder

Actions:
Copies data 3d INV2, UNI, & CEST acquisition dicoms to a new data directory and renames folders

Outputs: 
- Data is copied to a <subid_sesid> directory

#### 2a. convert2nifti.py
Hard-coded Inputs: paths 

#### 2b. convert2nifti.sh
Hard-coded Inputs: 
- Which files to convert; line 13-15 (for type in INV2 UNI cest b1map wassr none; do)

Actions: 
- Converts to nifti using:
    project/bbl_projects/apps/melliott/scripts/dicom2nifti.sh -u -F \

Outputs: 
- Niftis  

#### 3a. Preprocess_cest.sh
Inputs: path to niftis, output path
Actions: Calls process_GluCEST.py

#### 3b. Process_GluCEST.py
Inputs: 



