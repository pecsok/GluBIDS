import os

project = '/project/bbl_roalf_pecsokphd/projects/tms-glucest/data/cest'
scripts = '/project/bbl_roalf_pecsokphd/projects/glucest-rsfmri/github'
rawdata = os.path.join(project, 'rawdata')
out = os.path.join(project, 'rawdata')
convert = os.path.join(scripts, 'scripts_pyGluCEST', 'convert2nifti.sh')

for case in os.listdir(rawdata):
    cmd = ['bsub', convert, rawdata, out, case]
    os.system(' '.join(cmd))

