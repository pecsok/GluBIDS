#!/bin/sh

rawdata=$1
out=$2
case=$3

# make output directoreis for subject/session
if [ ! -d $out/$case ]; then
  mkdir $out/$case -p
fi


#if [ ! -e $out/$case/$case-INV2.nii ]; then
  # convert and bias correct INV2 and UNI
 for type in INV2 UNI cest b1map wassr none; do
     # Unzip dicoms
     #unzip $rawdata/"$case"/S*$type*/*dicom.zip
     /project/bbl_projects/apps/melliott/scripts/dicom2nifti.sh -u -F \
     $out/"$case"/$case-$type.nii \
     $rawdata/"$case"/S*$type*/*dcm
 done

#  for type in INV2 UNI cest none; do
  #  newnifti=$(dcm2niix -o $out/$case/ -v 0 -b y $rawdata/$case/S*_mark_*$type*/*dcm)

   # newnifti_name='$case-$type'
    #mv $newnifti $newnifti_name.nii
  #done

#fi
