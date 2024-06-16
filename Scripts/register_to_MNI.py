import os
import sys
import argparse


def register_to_MNI(data, out, structural, atlas, resolution):

    # load c function calc_b0_map and define return type and argument types
    script_path = os.path.realpath(__file__)
    script_path = os.path.dirname(script_path)  # 'scripts' folder
    script_path = os.path.join(script_path, 'register_to_MNI.sh')

    # iterate through the rawdata path and submit a job for each session
    for case in os.listdir(data):
        case_dir = os.path.join(data, case)
        
        #if os.path.isdir(os.path.join(case_dir, 'MNI_transforms')):
            #continue
        
        cmd = ['bsub', script_path, case, data, out, 
               structural, atlas, resolution]
        os.system(' '.join(cmd))
    
    return


def main():
    # Set up the argparser
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')

    # Add the argument for the input directory
    required.add_argument(
        "-d", "--data",
        help="Path to subject/session folders with brain masked INV2 and UNI",
        metavar="",
        required=True)
    
    required.add_argument(
        "-o", "--out",
        help="Path to output directory",
        metavar="",
        required=True)

    required.add_argument(
        "-s", "--structural",
        help="Type of structural image (ex. UNI, INV2, or mprage)",
        metavar="",
        required=True)

    required.add_argument(
        "-a", "--atlas",
        help="Path to atlases",
        metavar="",
        required=True)

    required.add_argument(
        "-r", "--resolution",
        help="Pixel dimension of the structural image",
        metavar="",
        required=True)


    # Parse the arguments
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    # run hdbet_seg
    register_to_MNI(args.data, args.out, args.structural,
                    args.atlas, args.resolution)


if __name__ == "__main__":
    main()

