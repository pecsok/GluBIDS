import os
import sys
import argparse
import warnings
from re import search
from write_images import write_image
from read_images import load_all_images
from calc_b0b1_map import calc_b0_map, calc_b1_map
from create_cest_masks import create_masks
from calculate_cest import calc_raw_cest, calc_b0_cest, calc_b0b1_cest


# helper function for main() that checks all the paths given by user
def check_paths(dir_list):

    for d in dir_list:
        if d is not None and not os.path.isdir(d):
            raise ValueError('Path {d} does not exist.'.format(d=repr(d)))

    return


# get the dataset's mask directory
def get_mask_dir(use_matlab_masks, use_mni_mask, mask, seg_mask_out_dir,
                 indir, outdir):

    # set mask_type to determine what kind of mask we are creating
    # if it's not provided
    if use_matlab_masks:
        mask_type = 'mat'
    elif use_mni_mask:
        mask_type = 'mni'
    # not using mat masks and using other user provided masks in args.mask
    elif len([x for x in os.listdir(mask) if 'mask.nii' in x]) == 1:
        mask_type = 'exists'
    # need to create mask (hdbet_seg) if not using mat, mni, nor existing masks
    else:
        mask_type = 'hdbet'

    # get the dataset's mask directory
    # using user provided masks that is NOT mat masks
    if mask_type == 'exists' or mask_type == 'mat':
        mask_dir = mask

    else:
        mask_dir = seg_mask_out_dir

    # find segmentation directory
    if len([x for x in os.listdir(mask) if 'seg.nii' in x]) == 1:
        seg_type = 'exists'
        seg_dir = mask
    else:
        seg_type = 'create'
        seg_dir = seg_mask_out_dir

    return [mask_type, mask_dir, seg_type, seg_dir]


# helper function to write out the estimated b0, b1, and b0b1 corrected CEST.
def write_all(b0map, b1map, b0b1_cest_map, sub_out,
              subject, mode, affine_dtype):

    # use the reference image data type and affine matrix
    # to write the proecessed images
    [affine, dtype] = affine_dtype
    write_image(b0b1_cest_map,
                os.path.join(sub_out, subject + '-' + 'B0B1CESTMAP.nii'),
                affine, dtype, mode)

    write_image(b0map,
                os.path.join(sub_out, subject + '-' + 'B0MAP.nii'),
                affine, dtype, mode)

    write_image(b1map,
                os.path.join(sub_out, subject + '-' + 'B1MAP.nii'),
                affine, dtype, mode)

    return


def process_GluCEST(data_dir, out_dir, mask_dir, mask_type, seg_type,
                    mode, dim, CESTppm, case=None):
    """ Process GluCEST and write resulting images


    Parameters:

    data_dir (str): path to directory that stores 2d CEST, B0, B1, 2d T1 maps by subject.
    out_dir (str): path to directory that write output images to.
    mode (str): "hippocampus" or "syrp" or "3d".
    dim (int): dimension of GluCEST image: 2 or 3.
    CESTppm (int): CEST ppm.
    mask_dir (str): Path where all subject's brain masks and segmentation maps are saved.
    use_matlab_masks (bool): True or False specifying if .mat file should be \
                             used for brain masks and segmentation maps. Assumes \
                             that the .mat file is in the data directory (data_dir).


    """

    # iterate through each subject in data_dir
    for subject in os.listdir(data_dir):

        if not search(".+_.+[^hippo]", subject) or \
                (case is not None and subject != case):
            continue

        # make sure output path exists
        sub_out = os.path.join(out_dir, subject)
        if not os.path.isdir(sub_out):
            os.mkdir(sub_out)

        # get the absolute path of subject and iterate to parse each file type
        sub_data = os.path.join(data_dir, subject)

        # get the absolute path of subject's mask directory
        sub_masks = os.path.join(mask_dir, subject)

        # create subject's mask and seg map if they don't exist
        create_masks(sub_data, sub_out, sub_masks, mask_type, seg_type)

        # make dictionary of dicom folder names for 2d and 3d
        # this is used to pass into load_all_images to parse sub_data directory
        fdict = {}
        if dim == 2:
            fdict.update({"cest_dcm": "moco_cest"})
            fdict.update({"wassr_dcm": "moco_wassr"})
            fdict.update({"b1_dcm": "moco_b1map"})
            fdict.update({"t1_dcm": None})
            fdict.update({"none_nifti": "none.nii"})

        else:
            fdict.update({"cest_dcm": "prep_tfl_3d_cest"})
            fdict.update({"wassr_dcm": "prep_tfl_3d_wassr"})
            fdict.update({"b1_dcm": "prep_tfl_3d_b1map"})
            fdict.update({"t1_dcm": "t1_images"})
            fdict.update({"none_nifti": "sat.nii"})

        # load all the images and return error if it fails
        # TODO: edit the downstream code for load_all_images to include maskZ_type
        #import pdb; pdb.set_trace()
        [cest_info, t1_map, wassr_info, b1_info,
         mask, seg, affine_dtype] = load_all_images(sub_data, sub_masks,
                                                    mask_type, mode,
                                                    sub_out, fdict, dim)

        # calculate B0 map
        b0map = calc_b0_map(wassr_info["WASSR_ppm_list"],
                            wassr_info["WASSR_pos_img"],
                            wassr_info["WASSR_neg_img"],
                            mask, dim, 0.0050)

        # calculate B1 map
        b1map = calc_b1_map(b1_info["B1_img1"], b1_info["B1_img2"],
                            b1_info["B1_img3"], mask, b1_info["alpha"],
                            dim, b1_info["B1_hdr"]["reps"])

        # calclate B0B1 corrected CEST map
        if dim == 2:
            # TODO:
            # calculate raw cest map for 2d (but it is not written nor used..?)
            raw_cest = calc_raw_cest(CESTppm,
                                     cest_info["CEST_ppm_list"],
                                     cest_info["CEST_neg_images"],
                                     cest_info["CEST_pos_images"],
                                     mask)
            # calculate b0 corrected cest
            b0_cest_map = calc_b0_cest(CESTppm, cest_info["CEST_ppm_list"],
                                       cest_info["CEST_pos_images"],
                                       cest_info["CEST_neg_images"],
                                       b0map, mask, dim)
            # calculate B0B1 corrected CEST maps
            b0b1_cest_map = calc_b0b1_cest(b1map, mask, seg,
                                           mode, b0_cest_map)

        else:
            # calculate b0 corrected cest
            [b0_corr_pos, b0_corr_neg] = calc_b0_cest(CESTppm,
                                                      cest_info["CEST_ppm_list"],
                                                      cest_info["CEST_pos_images"],
                                                      cest_info["CEST_neg_images"],
                                                      b0map, mask, dim)
            # calculate B0B1 corrected CEST maps
            t1_map = t1_map * mask
            b0b1_cest_map = calc_b0b1_cest(b1map, mask,
                                           cest_b1=cest_info["CEST_b1"],
                                           b0_corr_pos=b0_corr_pos,
                                           b0_corr_neg=b0_corr_neg,
                                           t1_map=t1_map)

        write_all(b0map, b1map, b0b1_cest_map,
                  sub_out, subject, mode, affine_dtype)

    return


def main():
    # Set up the argparser
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Add the argument for the input directory
    required.add_argument(
        "-i", "--input",
        help="Path to input directory with subject folders",
        metavar="",
        required=True)
    required.add_argument(
        "-o", "--output",
        help="Path to output directory",
        metavar="",
        required=True)
    required.add_argument(
        "-t", "--type",
        help="Type of data (syrp, hippocampus, or 3d)",
        metavar="",
        required=True)
    required.add_argument(
        "-d", "--dim",
        help="Dimension of data (2 or 3)",
        metavar="",
        required=True)
    optional.add_argument(
        "-m", "--mask",
        help="Directory that contains brain masks and segmentation maps. \n" +
             "Default is input directory and either '-b' or '-s' must be set " +
             " to 'exist' to use existing maps.",
        metavar="")
    optional.add_argument(
        "-p", "--cestppm",
        metavar="",
        help="CESTppm: default is 3")
    # optional.add_argument(
    #     "-s", "--seg_mask_out_dir",
    #     metavar="",
    #     help="Output direcotry for 3D and 2D brain masks and segmentation maps." +
    #          " Default is --output.")
    optional.add_argument(
        "-b", "--mask_type",
        metavar="",
        help="Select what kind of brain mask to use: \n" +
             "'exist': brain mask is in --mask \n" +
             "'hdbet': use built-in HDBET to brain mask \n" +
             "'mni': use MNI brain mask \n " +
             "'mat': use .mat file for brain mask \n" +
             "Default is 'hdbet'")
    optional.add_argument(
        "-s", "--seg_type",
        metavar="",
        help="Select what kind of brain mask to use: \n" +
             "'exist': segmentation map is in --mask \n" +
             "'create': use built-in FSL FAST+FIRST for segmentation \n" +
             "Default is 'create'.")
    optional.add_argument(
        "-c", "--case",
        metavar="",
        help="Case label (ex. <subject>_<session> in data directory)")
    parser._action_groups.append(optional)

    # Parse the arguments
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    # check the values of the arguments
    if args.type not in ['syrp', 'hippocampus', '3d'] or \
       args.dim not in ['2', '3'] or \
       (args.mask_type is not None and
            args.mask_type not in ['exist', 'hdbet', 'mni', 'mat']) or \
       (args.seg_type is not None and
            args.seg_type not in ['exist', 'create']):

        parser.print_help()
        sys.exit(2)

    # set variables
    args.dim = int(args.dim)

    # if user didn't provide args.mask, default is the input directory
    if args.mask is None:
        args.mask = args.input

    if args.cestppm:
        args.cestppm = float(args.cestppm)
    else:
        args.cestppm = 3.0

    if args.mask_type is None:
        args.mask_type = 'hdbet'

    if args.seg_type is None:
        args.seg_type = 'create'

    if args.mask is None and \
            (args.seg_type == 'exist' or args.mask_type == 'exist'):
        args.mask = args.input

    # check if paths exist
    check_paths([args.input, args.output, args.mask])

    # process GluCEST
    process_GluCEST(args.input, args.output, args.mask,
                    args.mask_type, args.seg_type,
                    args.type, args.dim, args.cestppm, args.case)

    return


if __name__ == "__main__":
    main()
