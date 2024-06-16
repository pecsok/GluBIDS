import os
from re import search

project = '/project/bbl_roalf_longglucest'
data = os.path.join(project, 'TMS_analysis/data')
dicoms = os.path.join(project, 'TMS_analysis/rawdata')

for sub in os.listdir(dicoms):

    if not search(r'^\d+', sub):
        continue

    sub_dicom = os.path.join(dicoms, sub)

    for ses in os.listdir(sub_dicom):

        if not search(r'^\d+', ses):
            continue

        ses_dicom = os.path.join(sub_dicom, ses)
        #ses_data = os.path.join(data, sub + '_' + ses)

        ses_data_onsite = os.path.join(data, sub + '_' + ses, 'onsite')
        ses_data_offsite = os.path.join(data, sub + '_' + ses, 'offsite')

        if not os.path.isdir(ses_data_onsite) and not os.path.isdir(ses_data_offsite):
            os.makedirs(ses_data_onsite)
            os.makedirs(ses_data_offsite)

            # make case folder in data
            #required_files = ['*mp2rage_mark_ipat3_0.80mm_INV2',
                              #'*mp2rage_mark_ipat3_0.80mm_UNI_Images',
                              #'*prep_moco_b1map*', '*prep_moco_cest*',
                              #'*prep_moco_none*',
                              #'*prep_moco_wassr*']

            required_files = ['*mp2rage_sag_3D_1x1x2mm_4nex_nsIR_250V_INV2',
                              '*mp2rage_sag_3D_1x1x2mm_4nex_nsIR_250V_T1_Images',
                              '*mp2rage_sag_3D_1x1x2mm_4nex_nsIR_250V_UNI_Images',
                              '*prep_tfl_3D_b1map_sag_usethis', '*prep_tfl_3D_cest_sag_usethis_210Hz',
                              '*prep_tfl_3D_sat_sag_usethis',
                              '*prep_tfl_3D_wassr_sag_usethis']

            for f in required_files:
                cmd = ['cp', '-r',
                       os.path.join(ses_dicom, 'onsite', f),
                       ses_data_onsite]
                cmd2 = ['cp', '-r', os.path.join(ses_dicom, 'offsite', f), ses_data_offsite]

                os.system("printf 'copying " + sub + '_' + ses + "\n'")
                os.system(' '.join(cmd))
                os.system(' '.join(cmd2))
