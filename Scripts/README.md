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

*****
Phase 2: Preprocessing
*****

#### 3a. Preprocess_cest.sh
Inputs: path to niftis, output path
Actions: Calls process_GluCEST.py

#### 3b. Process_GluCEST.py
Inputs: 

Actions: 
-Applies segmentation mask.
-Calculates B0 map (line 148-151) 
   - Library contains calc_b0_map.c, *.h, and *.so files
   - Input must be dicoms
-Calculates B1 map
-Calculates corrected B0B1 maps

Outputs: ***

#### 4a. register_to_MNI.py
Inputs: 
-Paths (input, output, atlases)
-Anatomical scan used for segmentation (UNI or INV2)
-Resolution (our data is 0.82)

Actions: Calls register_to_MNI.sh

#### 4b. register_to_MNI.sh
Inputs:
see above.

Actions:
-Checks to see if B0 and B1 correction and masking occurred successfully
-If not, uses old tools to perform (e.g., N4BiasCorrect; fslmaths -mas; antsRegistrationSyN.sh)

Outputs after steps 3 and 4:

*****
Phase 3: Data extraction
*****

#### 5a. process_cest.py
Inputs: 
paths (preprocessed dir, postprocessed dir, atlases, logs)
Segmentation input (UNI)
Resolution (0.82)

Actions: 
Calls process_cest.sh

#### 5b. process_cest.sh
Inputs: 
see above.

Actions: 
-Mask cest data and b0 maps
-Apply atlases

#### 6. extract_rois.sh
Inputs: 
(preprocessed dir, postprocessed dir, output dir)

Actions:
Loops through ROIs of designated atlas and extracts data into tsvs.
