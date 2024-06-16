#!/bin/bash

#This script calculates GluCEST contrast and gray matter density measures

#######################################################################################################
## DEFINE PATHS ##

cest=$1 #path to processed GluCEST data
data=$2
outputpath=$3

# HARVARD OXFORD

if ! [ -d $outputpath/INV2 ]; then
  mkdir $outputpath/INV2
fi

if ! [ -d $outputpath/UNI ]; then
  mkdir $outputpath/UNI
fi

for str in INV2 UNI; do

  for atlas in cort sub; do
    touch $outputpath/$str/all-GluCEST-HarvardOxford-$atlas-measures_$str.tsv
    echo "Subject	HarvardOxford_${atlas}_CEST_mean	HarvardOxford_CEST_numvoxels	HarvardOxford_CEST_SD" >> \
       	 $outputpath/$str/all-GluCEST-HarvardOxford-$atlas-measures_$str.tsv
  done

  for i in $(ls $cest); do
    case=${i##*/}
    echo "CASE: $case"
    mkdir $outputpath/$str/$case

    for atlas in cort sub; do
      # quantify GluCEST contrast for each participant
      3dROIstats -mask $cest/$case/atlases/$str/$case-cest-HarvardOxford-$atlas.nii.gz \
        -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
        $cest/$case/$case-GluCEST.nii.gz >> $outputpath/$str/$case/$case-HarvardOxfordROI-GluCEST-$atlas-measures_$str.tsv
      #format participant-specific csv
      sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-HarvardOxfordROI-GluCEST-$atlas-measures_$str.tsv
      cut -f2-3 --complement $outputpath/$str/$case/$case-HarvardOxfordROI-GluCEST-$atlas-measures_$str.tsv >> \
        $outputpath/$str/$case/tmp.tsv
      mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-HarvardOxfordROI-GluCEST-$atlas-measures_$str.tsv

      # quantify GluCEST contrast for each participant (whole subcortical and cortical)
      3dROIstats -mask $cest/$case/atlases/$str/$case-cest-HarvardOxford-$atlas-bin.nii.gz \
        -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
	$cest/$case/$case-GluCEST.nii.gz >> $outputpath/$str/$case/$case-HarvardOxford-GluCEST-$atlas-measures_$str.tsv
      #format participant-specific csv
      sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-HarvardOxford-GluCEST-$atlas-measures_$str.tsv
      cut -f2-3 --complement $outputpath/$str/$case/$case-HarvardOxford-GluCEST-$atlas-measures_$str.tsv >> \
        $outputpath/$str/$case/tmp.tsv
      mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-HarvardOxford-GluCEST-$atlas-measures_$str.tsv
      #enter participant GluCEST contrast data into master spreadsheet
      sed -n "2p" $outputpath/$str/$case/$case-HarvardOxford-GluCEST-$atlas-measures_$str.tsv >> \
        $outputpath/$str/all-GluCEST-HarvardOxford-$atlas-measures_$str.tsv

      # quantify 3d contrast for each participant
      3dROIstats -mask $cest/$case/atlases/$str/$case-HarvardOxford-$atlas.nii.gz \
	-zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat \
        $data/$case/$case-INV2_corrected.nii.gz >> $outputpath/$str/$case/$case-HarvardOxfordROI-INV2-$atlas-measures_$str.tsv
      #format participant-specific csv
      sed -i 's/name/Subject/g' $outputpath/$str/$case/$case-HarvardOxfordROI-INV2-$atlas-measures_$str.tsv
      cut -f2-3 --complement $outputpath/$str/$case/$case-HarvardOxfordROI-INV2-$atlas-measures_$str.tsv >> \
        $outputpath/$str/$case/tmp.tsv
      mv $outputpath/$str/$case/tmp.tsv $outputpath/$str/$case/$case-HarvardOxfordROI-INV2-$atlas-measures_$str.tsv

    done
  done
done


