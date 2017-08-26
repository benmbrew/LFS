#!/bin/bash

home=/hpf/largeprojects/agoldenb/ben/Projects
project=${home}/LFS
test=${project}/Scripts/predict_age

# Run the jobs
for i in {1..10}; do # seed
 
echo "${test}/model_pipeline_funnorm_no_transform_full.R $i" | qsub -N "${i}" -l nodes=1:ppn=12,gres=localhd:1,vmem=20G,mem=20G,walltime=7:00:00:00 -o ${test}/Output -e ${test}/Error
 sleep 0.1
 done
   
