import os

project = '/project/bbl_roalf_longglucest'
scripts = '/project/bbl_roalf_longglucest/sandbox/ally'
rawdata = os.path.join(project, 'TMS_analysis/whole_brain_data')
out = os.path.join(project, 'TMS_analysis/whole_brain_data')
convert = os.path.join(scripts, 'scripts_copy', 'convert2nifti.sh')

for case in os.listdir(rawdata):
    cmd = ['bsub', convert, rawdata, out, case, 'onsite']
    cmd2 = ['bsub', convert, rawdata, out, case, 'offsite']
    os.system(' '.join(cmd))
    os.system(' '.join(cmd2))

