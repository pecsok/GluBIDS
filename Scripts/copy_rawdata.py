import os
from re import search

project = '/project/bbl_roalf_pecsokphd/projects/tms-glucest/data/cest'
data = os.path.join(project, 'rawdata')
dicoms = os.path.join(project, '7T_data')

for sub in os.listdir(dicoms):

    if not search(r'^\d+', sub):
        continue

    sub_dicom = os.path.join(dicoms, sub)

    for ses in os.listdir(sub_dicom):

        if not search(r'^\d+', ses):
            continue

        ses_dicom = os.path.join(sub_dicom, ses)
        ses_data = os.path.join(data, sub + '_' + ses)

        if not os.path.isdir(ses_data):
            os.mkdir(ses_data)

            # make case folder in data
            """
            # Required files for 2D:
            required_files = ['*mp2rage_mark_ipat3_0.80mm_INV2',
                              '*mp2rage_mark_ipat3_0.80mm_UNI_Images',
                              '*prep_moco_b1map*', '*prep_moco_cest*',
                              '*prep_moco_none*',
                              '*prep_moco_wassr*']
            """

            required_files = ['*mp2rage_mark_ipat3_0.80mm_INV2',
                              '*mp2rage_3D_1x1x2mm_4nex_nsIR_250V_T1_Images',
                              '*mp2rage_mark_ipat3_0.80mm_UNI_Images',
                              '*prep_tfl_3D_b1map_hippo_usethis', '*prep_tfl_3D_cest_hippo_usethis_210Hz',
                              '*prep_tfl_3D_sat_hippo_usethis',
                              '*prep_tfl_3D_wassr_hippo_usethis']

            for f in required_files:
                cmd = ['cp', '-r',
                       os.path.join(ses_dicom, f),
                       ses_data]

                os.system("printf 'copying " + sub + '_' + ses + "\n'")
                os.system(' '.join(cmd))
