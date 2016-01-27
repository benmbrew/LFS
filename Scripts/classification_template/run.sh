#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben/Projects
project=${home}/LFS
test=${project}/Scripts/classification_template

# Run the jobs
  echo "Rscript ${test}/main.R $i" | qsub -N "${i}" -l nodes=2:ppn=35,gres=localhd:1,vmem=20G,mem=20G,walltime=4:00:00:00 -o ${project}/Output -e ${project}/Error

