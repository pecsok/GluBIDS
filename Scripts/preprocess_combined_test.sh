#!/bin/bash

scripts_path=$1

python $scripts_path/copy_rawdata.py &

wait

echo "finished copy_rawdata.py"

python $scripts_path/convert2nifti.py &
wait
echo "finished convert2nifti.py"

./$scripts_path/preprocess_cest.sh /path/2/data /path/2/output &
wait
echo "finished preprocess_cest.sh"

