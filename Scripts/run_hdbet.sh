#!/bin/sh

case=$1
case_rawdata=$2
case_out=$3
hdbet_image=$4
struc_img=$5

# convert and bias correct INV2
if [ ! -e $case_out/${struc_img}_corrected.nii.gz ]; then
  N4BiasFieldCorrection -d 3 -i $case_rawdata/${struc_img}.nii* \
    -o $case_out/${struc_img}_corrected.nii.gz
fi

# brain mask the INV2 using hdbet
if [ ! -e $case_out/${struc_img}-hdbet_mask.nii.gz ]; then

  singularity run $hdbet_image \
    -i $case_out/${struc_img}_corrected.nii.gz \
    -o $case_out/${struc_img}-hdbet.nii.gz \
    -device cpu \
    -mode fast \
    -tta 0

fi

printf "Finished run_hdbet.sh\n"
