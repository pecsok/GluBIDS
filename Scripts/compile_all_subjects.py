import os
import sys
import argparse
import pandas as pd
from xml.dom import minidom


def compile_rois(data, out, xml, group, str_type):

    # make an empty dataframe that will have all subjects' roi values
    group = pd.read_csv(group)
    group.columns = ['Subject', 'group']
    for roi in ["cort", "sub"]:
        if roi == 'cort':
            idx = os.path.join(xml, 'HarvardOxford-Cortical.xml')
        else:
            idx = os.path.join(xml, 'HarvardOxford-Subcortical.xml')

        # read the xml file that has all the labels for each JHU roi
        idx = minidom.parse(idx)
        labels = idx.getElementsByTagName('label')

        for img in [str_type, 'GluCEST']:
            # read in the case's ROI values for ROI roi
            roi_table = img + '-HarvardOxfordROI-' + roi + '-Measures_' + str_type + '.tsv'
            roi_table = pd.read_table(os.path.join(out, str_type, roi_table), delimiter='\t')

            # get colnames of roi_table
            roi_columns = list(roi_table.columns)
            new_columns = ['Subject']
            for col in roi_columns:
                if col == 'Subject':
                    continue
                roi_parts = col.split('_')
                roi_measure = roi_parts[0]
                roi_num = roi_parts[-1]
                roi_num = int(roi_num) - 1

                # get the label
                roi_label = labels[roi_num].firstChild.data
                new_columns.append(roi_label + ' ' + roi_measure)

            # get subject labels
            subject = roi_table["Subject"]
            new_subject = []
            for sub in subject:
                sub = sub.split('/')[-1]
                sub = sub.split('-')[0]
                new_subject.append(sub)

            roi_table['Subject'] = new_subject
            roi_table.columns = new_columns
            roi_table = pd.merge(group, roi_table, how = 'right')
            roi_table.to_csv(os.path.join(out, str_type, 
                                          'all_subs_' + img + '_' + roi + '_rois_' + str_type + '.csv'))

    return


def main():
    # Set up the argparser
    parser = argparse.ArgumentParser()
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
        "-x", "--xml",
        help="path to HarvardOxford xml files",
        metavar="",
        required=True)
    required.add_argument(
        "-g", "--group",
        help="path to CEST_Groups.csv file",
        metavar="",
        required=True)

    # Parse the arguments
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)


    compile_rois(args.input, args.output, args.xml, args.group, 'UNI')
    #compile_rois(args.input, args.output, args.xml, args.group, 'INV2')

    return


if __name__ == "__main__":
    main()
