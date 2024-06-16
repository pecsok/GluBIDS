#!/bin/bash

#This script calculates GluCEST contrast and gray matter density measures

#######################################################################################################
## DEFINE PATHS ##

cest=$1 #path to processed GluCEST data
data=$2
outputpath=$3

# SCHAEFER 2018

if ! [ -d $outputpath/INV2 ]; then
  mkdir $outputpath/INV2
fi

if ! [ -d $outputpath/UNI ]; then
  mkdir $outputpath/UNI
fi

for str in INV2 UNI; do
  touch $outputpath/$str/all-GluCEST-Schaefer2018-100P-17N-measures_$str.tsv
  echo "Subject	Schaefer2018_CEST_mean	Schaefer2018_CEST_numvoxels	Schaefer2018_CEST_SD" >> \
       $outputpath/$str/all-GluCEST-Schaefer2018-100P-17N-measures_$str.tsv

  for i in $(ls $cest); do
    case=${i##*/}
    echo "CASE: $case"
    mkdir $outputpath/$str/$case

    # quantify GluCEST contrast for each participant
    3dROIstats -mask $cest/$case/atlases/$str/$case-cest-Schaefer2018-100P-17N.nii.gz \
      -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
      $cest/$case/$case-GluCEST.nii.gz >> $outputpath/$str/$case/$case-Schaefer2018ROI-GluCEST-100P-17N-measures_$str.tsv
    # format participant-specific csv
    sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-Schaefer2018ROI-GluCEST-100P-17N-measures_$str.tsv
    cut -f2-3 --complement $outputpath/$str/$case/$case-Schaefer2018ROI-GluCEST-100P-17N-measures_$str.tsv >> \
      $outputpath/$str/$case/tmp.tsv
    mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-Schaefer2018ROI-GluCEST-100P-17N-measures_$str.tsv

    3dROIstats -mask $cest/$case/atlases/$str/$case-cest-Schaefer2018-100P-17N-bin.nii.gz \
      -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
      $cest/$case/$case-GluCEST.nii.gz >> $outputpath/$str/$case/$case-Schaefer2018-GluCEST-100P-17N-measures_$str.tsv
    # format participant specific csv
    sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-Schaefer2018-GluCEST-100P-17N-measures_$str.tsv
    cut -f2-3 --complement $outputpath/$str/$case/$case-Schaefer2018-GluCEST-100P-17N-measures_$str.tsv >> \
      $outputpath/$str/$case/tmp.tsv
    mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-Schaefer2018-GluCEST-100P-17N-measures_$str.tsv
    # enter participant GluCEST contrast data into master spreadsheet
    sed -n "2p" $outputpath/$str/$case/$case-Schaefer2018-GluCEST-100P-17N-measures_$str.tsv >> \
      $outputpath/$str/all-GluCEST-Schaefer2018-100P-17N-measures_$str.tsv

    # quantify 3d contrast for each participant
    3dROIstats -mask $cest/$case/atlases/$str/$case-Schaefer2018-100P-17N.nii.gz \
    -zerofill Nan -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
    $data/$case/$case-INV2_corrected.nii.gz >> $outputpath/$str/$case/$case-Schaefer2018ROI-INV2-100P-17N-measures_$str.tsv
    # format participant-specific csv
    sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-Schaefer2018ROI-INV2-100P-17N-measures_$str.tsv
    cut -f2-3 --complement $outputpath/$str/$case/$case-Schaefer2018ROI-INV2-100P-17N-measures_$str.tsv >> \
      $outputpath/$str/$case/tmp.tsv
    mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-Schaefer2018ROI-INV2-100P-17N-measures_$str.tsv

    for network in Cont Default DorsAttn Limbic SalVentAttn SomMot Vis
    do
      3dROIstats -mask $cest/$case/atlases/$str/$case-2d-s100_7-${network}.nii.gz \
        -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
        $cest/$case/$case-GluCEST.nii.gz >> $outputpath/$str/$case/$case-2dROI-GluCEST-s100_7-${network}-measures_$str.tsv
      #format participant-specific csv
      sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-2dROI-GluCEST-s100_7-${network}-measures_$str.tsv
      cut -f2-3 --complement $outputpath/$str/$case/$case-2dROI-GluCEST-s100_7-${network}-measures_$str.tsv >> \
        $outputpath/$str/$case/tmp.tsv
      mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-2dROI-GluCEST-s100_7-${network}-measures_$str.tsv

      3dROIstats -mask $cest/$case/atlases/$str/$case-2d-s100_7-${network}-bin.nii.gz \
        -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
        $cest/$case/$case-GluCEST.nii.gz >> $outputpath/$str/$case/$case-2d-GluCEST-s100_7-${network}-measures_$str.tsv
      # format participant specific csv
      sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-2d-GluCEST-s100_7-${network}-measures_$str.tsv
      cut -f2-3 --complement $outputpath/$str/$case/$case-2d-GluCEST-s100_7-${network}-measures_$str.tsv >> \
        $outputpath/$str/$case/tmp.tsv
      mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-2d-GluCEST-s100_7-${network}-measures_$str.tsv
      # enter participant GluCEST contrast data into master spreadsheet
      sed -n "2p" $outputpath/$str/$case/$case-2d-GluCEST-s100_7-${network}-measures_$str.tsv >> \
        $outputpath/$str/all-GluCEST-s100_7-${network}-measures_$str.tsv
    done
  done
done
