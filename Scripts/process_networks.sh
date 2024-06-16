#!/bin/bash

## DEFINE PATHS ##
structural=$1
pre=$2
post=$3
atlas=$4
case=$5
str=$6


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
      $post/$case/atlases/$str/${case}-2d-cest-s100_7-$network.nii

    gzip $post/$case/atlases/$str/${case}-2d-cest-s100_7-$network.nii

    # Multiply masked CEST file with binary tissue mask
    fslmaths $post/$case/atlases/$str/${case}-2d-cest-s100_7-$network.nii.gz \
      -mul $post/$case/$case-tissuemap-bin.nii.gz \
      $post/$case/atlases/$str/$case-2d-cest-s100_7-$network.nii.gz

    fslmaths $post/$case/atlases/$str/$case-2d-cest-s100_7-$network.nii.gz \
      -bin $post/$case/atlases/$str/$case-2d-cest-s100_7-${network}-bin.nii
done

#################################################################################################################
# Schaefer100 7 networks Method 2: Generate masks as you go
# Warp entire atlas to subject space
# antsApplyTransforms -d 3 -r $structural/$case/$case-UNI-masked.nii.gz \
# -i $structural/MNI_Templates/Schaefer100_7network_atlas-0.8mm.nii.gz \
# -n MultiLabel -o $structural/$case/atlases/${case}-s100_7_atlas.nii.gz \
# -t [$structural/$case/MNI_transforms/$case-UNIinMNI-0GenericAffine.mat,1] \
# -t $structural/$case/MNI_transforms/$case-UNIinMNI-1InverseWarp.nii.gz
# Extract2slice
# /project/bbl_projects/apps/melliott/scripts/extract_slice2.sh \
# -MultiLabel $structural/$case/atlases/${case}-s100_7_$network.nii.gz \
# $cest/$case/orig_data/$case-B0B1CESTMAP.nii \
# $cest/$case/atlases/${case}-2d-s100_7_atlas.nii
	
# gzip $cest/$case/atlases/${case}-2d-s100_7_atlas.nii
# Extract networks by thresholding atlas to include only the numerical values of the parcels we want.
# List of networks and thresholds
# networks = [Vis SomMot DorsAttn Limbic SalVentAttn Cont Default]
# Fyi: all_thresholds = [Vis: 1-9, 51-58; SomMot: 10-15, 59-66; DorsAttn: 16-23, 67-73; SalVentAttn: 24-30, 74-78; Limbic: 31-33, 79-80; Cont: 34-37, 81-89; Default: 38-50,90-100)
# Upper and lower thresholds for left hemisphere .. choose whichever is correct hemisphere for our slices!
# lth = [1 10 16 24 31 34 38]
# uth = [9 15 23 30 33 37 50]
# Upper and lower thresholds for right hemisphere 
# lth = [51 59 67 74 79 81 90]
# uth = [58 66 73 78 80 89 100]
# for i in 1:7
# do
# 	network = networks[i]
# 	upper = uth[i]
# 	lower = lth[i]
	# Multiply 2d Schaefer atlas file with tissue binary mask if not done previously
    # fslmaths $cest/$case/atlases/$case-2d-s100_7_atlas.nii.gz -mul $cest/$case/fast/$case-tissuemap-bin.nii.gz $cest/$case/atlases/$case-2d-s100_7-atlas.nii.gz
    # Threshold atlas to include only parcels in the network (need to figure out again if upper/lower bounds are inclusive
    # fslmaths $cest/$case/atlases/$case-2d-s100_7_atlas.nii.gz -thr ${lth} -uthr ${uth} $cest/$case/atlases/$case-2d-s100_7-${network}.nii.gz
    # Binarize mask
    # fslmaths $cest/$case/atlases/$case-2d-s100_7-${network}.nii.gz -bin $cest/$case/atlases/$case-2d-s100_7-${network}.nii.gz
    # Apply mask to CEST file and proceed with data extraction ...check this last step!
    # fslmaths $cest/$case/atlases/$case-2d-s100_7-atlas.nii.gz -mul $cest/$case/atlases/$case-2d-s100_7-${network} $cest/$case/atlases/$case-2d-RewardAtlas-Cortical.nii.gz
# done
