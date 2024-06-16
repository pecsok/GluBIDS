#!/bin/sh

data=$1
out=$2

for sub in $(ls $data)
do
  
  # if ! [ -d $out/$sub ]
  # then
      bsub python /project/bbl_roalf_longglucest/sandbox/ally/pyGluCEST_2.0/scripts/process_GluCEST.py \
      -i $data -o $out \
      -t syrp -d 2 \
      -m $out -b hdbet -s create \
      -c $sub
  # fi
done
