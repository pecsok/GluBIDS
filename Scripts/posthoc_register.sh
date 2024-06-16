#!/bin/bash

## DEFINE PATHS ##
structural=$1
pre=$2
post=$3
atlas=$4
log=$5
case=$6
method=$7
reso=$8
str=$9

#######################################################################################################
## IDENTIFY CASES FOR PROCESSING ##

echo "CASE: $case"


#check for structural data
if [ -e $structural/$case/MNI_transforms/$case-${str}inMNI-Warped.nii.gz ]
then
echo "Structural Data exists for $case"
sleep 1.5
else
echo "Oh No! Structural Data is missing. Cannot process CEST! Run register_to_MNI.sh first."
return
fi

#check for GluCEST GUI data
if ( [ -d $pre/$case/*WASSR_B0MAP2D ] && \
[ -d $pre/$case/*B1MAP2D ] && \
[ -d $pre/$case/*B0B1CESTMAP2D ] ) || \
( [ -e $pre/$case/$case-B0B1CESTMAP.nii ] && \
[ -e $pre/$case/$case-B0MAP.nii ] && \
[ -e $pre/$case/$case-B1MAP.nii ] )
then

echo "CEST GUI Data exists for $case"
sleep 1.5
else
echo "Oh No! CEST GUI Data is missing. Cannot process CEST! Analyze this case with CEST_2d_TERRA first."
sleep 1.5
fi

#if ! [ -d $post/$case ]
#then

logfile=$log/$case-cest.log
#{
echo "-------- Processing GluCEST data for $case ---------"

#######################################################################################################
## make directories and log files ##
#mkdir $post/$case -p
log_files=$post/$case/log_files #path to intermediate files. Remove for final script
#mkdir $log_files
mkdir $post/$case/atlases
#mkdir $post/$case/orig_data

#######################################################################################################
## CONVERT B0, B1, and B0B1-CORRECTED CEST FROM DCM TO NII ##
#######################################################################################################
## THRESHOLD B0 AND B1 MAPS ##

#threshold b0 from -1 to 1 ppm (relative to water resonance)
#fslmaths $pre/$case/$case-B0MAP.nii \
#  -add 10 \
#  $post/$case/$case-B0MAP-pos.nii.gz # make B0 map values positive to allow for thresholding with fslmaths
#######################################################################################################

## APPLY THRESHOLDED B0 MAP, B1 MAP, and TISSUE MAP (CSF removed) TO GLUCEST IMAGES ##

#exclude voxels with B0 offset greater than +- 1 pmm from GluCEST images
#######################################################################################################
## MASK THE PROCESSED GLUCEST IMAGE ##

#######################################################################################################
# clean up and organize, whistle while you work
#mv -f $post/$case/*masktmp* $log_files
#mv -f $post/$case/*.log $log_files
#mv -f $post/$case/$case-B0MAP-pos.nii.gz $log_files
#mv -f $post/$case/$case-B0MAP-thresh.nii.gz $log_files
#mv -f $post/$case/$case-B1MAP-thresh.nii.gz $log_files
# mv $post/$case/$case-B1MAP.nii $post/$case/orig_data
# mv $post/$case/$case-B0MAP.nii $post/$case/orig_data
# mv $post/$case/$case-B0B1CESTMAP.nii $post/$case/orig_data

#######################################################################################################
## REGISTER ATLASES TO UNI IMAGES AND GLUCEST IMAGES ##

#Harvard Oxford Atlasess
# for str in INV2 UNI; do
  
  if ! [ -d $post/$case/atlases/$str ]; then
    mkdir $post/$case/atlases/$str
  fi

  for roi in cort sub; do
    antsApplyTransforms -d 3 -r $structural/$case/${case}-${str}_masked.nii.gz \
      -i $atlas/HarvardOxford/HarvardOxford-$roi-maxprob-thr25-${reso}mm.nii.gz \
      -n MultiLabel \
      -o $post/$case/atlases/$str/${case}-HarvardOxford-$roi.nii.gz  \
      -t [$structural/$case/MNI_transforms/${case}-${str}inMNI-0GenericAffine.mat,1] \
      -t $structural/$case/MNI_transforms/${case}-${str}inMNI-1InverseWarp.nii.gz

    /project/bbl_projects/apps/melliott/scripts/extract_slice2.sh \
      -MultiLabel $post/$case/atlases/$str/${case}-HarvardOxford-$roi.nii.gz \
      $pre/$case/$case-B0B1CESTMAP.nii \
      $post/$case/atlases/$str/${case}-cest-HarvardOxford-$roi.nii

    gzip $post/$case/atlases/$str/${case}-cest-HarvardOxford-$roi.nii

    fslmaths $post/$case/atlases/$str/$case-cest-HarvardOxford-$roi.nii.gz \
      -mul $post/$case/$case-tissuemap-bin.nii.gz \
      $post/$case/atlases/$str/$case-cest-HarvardOxford-$roi.nii.gz

    fslmaths $post/$case/atlases/$str/$case-cest-HarvardOxford-$roi.nii.gz \
      -bin $post/$case/atlases/$str/$case-cest-HarvardOxford-$roi-bin.nii
  done

# done

#######################################################################################################
echo -e "\n$case SUCCESFULLY PROCESSED\n\n\n"
#} | tee "$logfile"
#else
#echo "$case is either missing data or already processed. Will not process"
#sleep 1.5
#fi

