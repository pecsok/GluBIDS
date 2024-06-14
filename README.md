
**GluBIDS**: A Preprocessing Pipeline for GluCEST Data

*****
Acquisitions
*****

Preliminary Naming Structure
```
sub-${bblid}/
  ses-${session}/
    anat/
      sub-${bblid}_ses-${session}_inv-1_part-mag_MP2RAGE.nii.gz
      sub-${bblid}_ses-${session}_inv-1_part-phase_MP2RAGE.json
      sub-${bblid}_ses-${session}_inv-2_part-mag_MP2RAGE.nii.gz
      sub-${bblid}_ses-${session}_inv-2_part-phase_MP2RAGE.nii.gz
      sub-${bblid}_ses-${session}_UNIT1.nii.gz  # T1w images
      sub-${bblid}_ses-${session}_UNIT1.json  # T1w images
      sub-${bblid}_ses-${session}_T1map.nii.gz  # Quantitative T1 image
      sub-${bblid}_ses-${session}_T1map.json  # Quantitative T1 image
    fmap
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_fieldmap.nii.gz  # prep_moco_none
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_fieldmap.json  # prep_moco_none
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_RB1map.nii.gz  # prep_moco_b1map
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_RB1map.json  # prep_moco_b1map
    cest
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_cest.nii.gz  # prep_moco_cest
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_cest.json  # prep_moco_cest
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_wassr.nii.gz  # prep_moco_wassr
      sub-${bblid}_ses-${session}_acq-3D_chunk-01_wassr.json  # prep_moco_wassr
```
