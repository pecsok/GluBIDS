#!/bin/sh

rawdata=$1
out=$2
case=$3
site=$4

# make output directoreis for subject/session
if [ ! -d $out/$case ]; then
  mkdir $out/$case -p
fi


#if [ ! -e $out/$case/$case-INV2.nii ]; then
  # convert and bias correct INV2 and UNI
 for type in INV2 UNI cest sat b1map wassr T1; do
     /project/bbl_projects/apps/melliott/scripts/dicom2nifti.sh -u \
     $out/"$case"/$case-$type.nii \
     $rawdata/"$case"/S*$type*/*dcm
 done

#  for type in INV2 UNI cest none; do
  #  newnifti=$(dcm2niix -o $out/$case/ -v 0 -b y $rawdata/$case/S*_mark_*$type*/*dcm)

   # newnifti_name='$case-$type'
    #mv $newnifti $newnifti_name.nii
  #done

#fi
