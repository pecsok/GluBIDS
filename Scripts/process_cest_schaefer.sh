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

if ! [ -d $post/$case ]
then

logfile=$log/$case-cest.log
{
echo "-------- Processing GluCEST data for $case ---------"

#######################################################################################################
## make directories and log files ##
mkdir $post/$case -p
log_files=$post/$case/log_files #path to intermediate files. Remove for final script
mkdir $log_files
mkdir $post/$case/atlases
mkdir $post/$case/orig_data

#######################################################################################################
## CONVERT B0, B1, and B0B1-CORRECTED CEST FROM DCM TO NII ##
if [[ "$method" == "matlab" ]]
then
  for seq in B0MAP B1MAP B0B1CESTMAP
  do
    /project/bbl_projects/apps/melliott/scripts/dicom2nifti.sh -u -r Y -F \
      $pre/$case/$case-$seq.nii $pre/$case/S*${seq}2D/*dcm
  done

fi
#######################################################################################################
## THRESHOLD B0 AND B1 MAPS ##

#threshold b0 from -1 to 1 ppm (relative to water resonance)
fslmaths $pre/$case/$case-B0MAP.nii \
  -add 10 \
  $post/$case/$case-B0MAP-pos.nii.gz # make B0 map values positive to allow for thresholding with fslmaths
fslmaths $post/$case/$case-B0MAP-pos.nii.gz \
  -thr 9 \
  -uthr 11 \
  $post/$case/$case-B0MAP-thresh.nii.gz #threshold from -1(+10=9) to 1(+10=11)
fslmaths $post/$case/$case-B0MAP-thresh.nii.gz \
  -bin $post/$case/$case-b0.nii.gz #binarize thresholded B0 map

#threshold b1 from 0.3 to 1.3
fslmaths $pre/$case/$case-B1MAP.nii \
  -thr 0.3 \
  -uthr 1.3 $post/$case/$case-B1MAP-thresh.nii.gz #threshold from 0.3 to 1.3
fslmaths $post/$case/$case-B1MAP-thresh.nii.gz \
  -bin $post/$case/$case-b1.nii.gz #binarize thresholded B1 map
#######################################################################################################

## APPLY THRESHOLDED B0 MAP, B1 MAP, and TISSUE MAP (CSF removed) TO GLUCEST IMAGES ##

#exclude voxels with B0 offset greater than +- 1 pmm from GluCEST images
fslmaths $pre/$case/$case-B0B1CESTMAP.nii \
  -mul $post/$case/$case-b0.nii.gz \
  $post/$case/$case-CEST_b0thresh.nii.gz

#exclude voxels with B1 values outside the range of 0.3 to 1.3 from GluCEST images
fslmaths $post/$case/$case-CEST_b0thresh.nii.gz \
  -mul $post/$case/$case-b1.nii.gz \
  $post/$case/$case-CEST_b0b1thresh.nii.gz

#exclude CSF voxels from GluCEST images
fslmaths $structural/$case/${case}_cestseg.nii -thr 2 $post/$case/$case-tissuemap.nii.gz
fslmaths $post/$case/$case-tissuemap.nii.gz -bin $post/$case/$case-tissuemap-bin.nii.gz
fslmaths $post/$case/$case-CEST_b0b1thresh.nii.gz -mul $post/$case/$case-tissuemap-bin.nii.gz $post/$case/$case-CEST-finalthresh.nii.gz

#######################################################################################################
## MASK THE PROCESSED GLUCEST IMAGE ##

fslmaths $pre/$case/$case-B1MAP.nii \
  -bin $post/$case/CEST-masktmp.nii.gz
fslmaths $post/$case/CEST-masktmp.nii.gz \
  -ero -kernel sphere 1 \
  $post/$case/CEST-masktmp-er1.nii.gz
fslmaths $post/$case/CEST-masktmp-er1.nii.gz \
  -ero -kernel sphere 1 \
  $post/$case/CEST-masktmp-er2.nii.gz
fslmaths $post/$case/CEST-masktmp-er2.nii.gz \
  -ero -kernel sphere 1 \
  $post/$case/$case-CEST-mask.nii.gz
fslmaths $post/$case/$case-CEST-finalthresh.nii.gz \
  -mul $post/$case/$case-CEST-mask.nii.gz \
  $post/$case/$case-GluCEST.nii.gz #final processed GluCEST Image
#######################################################################################################
# clean up and organize, whistle while you work
mv -f $post/$case/*masktmp* $log_files
mv -f $post/$case/*.log $log_files
mv -f $post/$case/$case-B0MAP-pos.nii.gz $log_files
mv -f $post/$case/$case-B0MAP-thresh.nii.gz $log_files
mv -f $post/$case/$case-B1MAP-thresh.nii.gz $log_files
# mv $post/$case/$case-B1MAP.nii $post/$case/orig_data
# mv $post/$case/$case-B0MAP.nii $post/$case/orig_data
# mv $post/$case/$case-B0B1CESTMAP.nii $post/$case/orig_data

#######################################################################################################
## REGISTER ATLASES TO UNI IMAGES AND GLUCEST IMAGES ##

# SCHAEFER 2018 100 parcels atlas

if ! [ -d $post/$case/atlases/$str ]; then
  mkdir $post/$case/atlases/$str
fi

antsApplyTransforms -d 3 -r $structural/$case/${case}-${str}_masked.nii.gz \
  -i $atlas/Schaefer2018/Schaefer2018_100Parcels_17Networks_${reso}mm.nii.gz \
  -n MultiLabel \
  -o $post/$case/atlases/$str/${case}-Schaefer2018-100P-17N.nii.gz \
  -t [$structural/$case/MNI_transforms/${case}-${str}inMNI-0GenericAffine.mat,1] \
  -t $structural/$case/MNI_transforms/${case}-${str}inMNI-1InverseWarp.nii.gz

/project/bbl_projects/apps/melliott/scripts/extract_slice2.sh \
  -MultiLabel $post/$case/atlases/$str/${case}-Schaefer2018-100P-17N.nii.gz \
  $pre/$case/$case-B0B1CESTMAP.nii \
  $post/$case/atlases/$str/${case}-cest-Schaefer2018-100P-17N.nii

gzip $post/$case/atlases/$str/${case}-cest-Schaefer2018-100P-17N.nii

fslmaths $post/$case/atlases/$str/${case}-cest-Schaefer2018-100P-17N.nii.gz \
  -mul $post/$case/$case-tissuemap-bin.nii.gz \
  $post/$case/atlases/$str/$case-cest-Schaefer2018-100P-17N.nii.gz

fslmaths $post/$case/atlases/$str/$case-cest-Schaefer2018-100P-17N.nii.gz \
  -bin $post/$case/atlases/$str/$case-cest-Schaefer2018-100P-17N-bin.nii

# Schaefer100 7 networks if inputting pre-made masks
for network in Cont Default DorsAttn Limbic SalVentAttn SomMot Vis
do
  # Warp mask to subject space
  antsApplyTransforms -d 3 -r $structural/$case/${case}-${str}_masked.nii.gz \
    -i $atlas/Schaefer2018/s100_7_${network}-0.8mm.nii.gz \
    -n MultiLabel \
    -o $post/$case/atlases/$str/${case}-s100_7_$network.nii.gz \
    -t [$structural/$case/MNI_transforms/${case}-${str}inMNI-0GenericAffine.mat,1] \
    -t $structural/$case/MNI_transforms/${case}-${str}inMNI-1InverseWarp.nii.gz

  # Extract2slice
  /project/bbl_projects/apps/melliott/scripts/extract_slice2.sh \
    -MultiLabel $post/$case/atlases/$str/${case}-s100_7_$network.nii.gz \
    $pre/$case/$case-B0B1CESTMAP.nii \
    $post/$case/atlases/$str/${case}-2d-s100_7_$network.nii

  gzip $post/$case/atlases/$str/${case}-2d-s100_7_$network.nii

  # Multiply masked CEST file with binary tissue mask
  fslmaths $post/$case/atlases/$str/${case}-2d-s100_7_$network.nii.gz \
    -mul $post/$case/$case-tissuemap-bin.nii.gz \
    $post/$case/atlases/$str/$case-2d-s100_7-$network.nii.gz

  fslmaths $post/$case/atlases/$str/$case-2d-s100_7-$network.nii.gz \
    -bin $post/$case/atlases/$str/$case-2d-s100_7-${network}-bin.nii
done


#######################################################################################################
echo -e "\n$case SUCCESFULLY PROCESSED\n\n\n"
} | tee "$logfile"
else
echo "$case is either missing data or already processed. Will not process"
sleep 1.5
fi

