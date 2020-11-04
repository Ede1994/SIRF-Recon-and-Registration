#!/bin/sh
# !! change path for reading AND writing !!

# start number
N=0

# for loop for converting dcm to nii
for file in /home/eric/Dokumente/PythonProjects/SIRF/UKL_data/EPI/20191009/1/MR* ; do

    miconv -rf dcm -type float $file /home/eric/Dokumente/PythonProjects/SIRF/UKL_data/EPI/20191009/1/epi_$N.nii

    ((N++))

    rm $file

done