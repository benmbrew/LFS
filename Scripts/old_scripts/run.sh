#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben/Projects
project=${home}/LFS
test=${project}/Scripts/bumphunter

# Run the jobs
echo "Rscript ${test}/get_probe_subset.R" | qsub -N main.R -l vmem=20G,mem=20G -o ${test}/Output -e ${test}/Error

