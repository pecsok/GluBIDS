case=$1
data=$2
out=$3
str=$4
atlas=$5
reso=$6


#id=$(echo $case | cut -d_ -f1,2)

# check output directory exists
if [ ! -d $out/$case/MNI_transforms ]; then
  mkdir $out/$case/MNI_transforms
fi

# check bias corrected structural image exists
if [ ! -e $out/$case/$case-${str}_corrected.nii.gz ]; then
  N4BiasCorrect -i $data/$case/$case-${str}.nii.gz \
    -o $data/$case/$case-${str}_corrected.nii.gz
fi

# check that the bias corrected structural image is skull-stripped
if [ ! -e $data/$case/$case-${str}_masked.nii.gz ]; then
  fslmaths $data/$case/$case-${str}_corrected.nii.gz \
    -mas $out/$case/$case-${str}*_mask.nii.gz \
    $data/$case/$case-${str}_masked.nii.gz
fi

# register the bias corrected & skull-stripped structural image to MNI
if [ ! -e $out/$case/MNI_transforms/$case-${str}inMNI-Warped.nii.gz ]; then
  #register brain masked INV2 to upsampled MNI T1 template
  antsRegistrationSyN.sh -d 3 \
    -f $atlas/MNI/MNI152_T1_${reso}mm_brain.nii.gz \
    -m $data/$case/$case-${str}_masked.nii.gz \
    -o $out/$case/MNI_transforms/$case-${str}inMNI-
fi


