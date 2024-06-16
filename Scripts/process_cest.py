import os
import sys
import argparse


def postprocess(script, data, pre, post, atlas, log, method, resolution, str_type):
    
    for case in os.listdir(data):

#        if not os.path.isdir(os.path.join(post, case)):  
            cmd = ['bsub', script, data, pre, post, 
                   atlas, log, case, method, resolution, str_type]
            
            os.system(' '.join(cmd))
    
    return


def main():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')

    required.add_argument(
        '-d', '--data',
        help = 'Path to data directory',
        required = True)
    required.add_argument(
        '-p', '--pre',
        help = 'path to preprocessed cest data',
        required = True)
    required.add_argument(
        '-o', '--post',
        help = 'Path to postprocessed output cest',
        required = True)
    required.add_argument(
        '-a', '--atlas',
        help = 'Path to atlases',
        required = True)
    required.add_argument(
        '-l', '--logs',
        help = 'Path to logs',
        required = True)
    required.add_argument(
        '-s', '--script',
        help = 'Path to process_cest.sh',
        required = True)
    required.add_argument(
        '-m', '--method',
        help = 'matlab or python',
        required = True)
    required.add_argument(
        '-r', '--resolution',
        help = 'Voxel resolution (ex. 0.7)',
        required = True)
    required.add_argument(
        '-t', '--str_type',
        help = 'Type of structural image (ex. INV2, UNI, mprage)',
        required = True)

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    postprocess(args.script, args.data, args.pre, args.post, 
                args.atlas, args.logs, args.method, args.resolution, args.str_type)

    return


if __name__ == "__main__":
    main()


