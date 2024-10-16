#!/bin/bash

# slurm settings
#SBATCH --job-name=MBA
#SBATCH --mem=4G
#SBATCH --partition=luna-cpu-long
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=24
#SBATCH --time=02-0:00:00
#SBATCH --nice=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --output=MBA_%A.out


# running containerized version of Matrix-Based Analysis Program through Bayesian Multilevel Modeling available in AFNI
#https://afni.nimh.nih.gov/pub/dist/doc/program_help/MBA.html
# Chen et al. 2019 paper - https://doi.org/10.1101/459545

export APPTAINER_BIND="/data/anw/anw-work,/scratch"
workdir=/data/anw/anw-work/NP/projects/data_chris/CORE/stats/MBA

for modality in dwi func; do 


# -  percentage improvement
modelID=${modality}_perc_improv

mkdir -p ${workdir}/${modelID}
cd ${workdir}/${modelID}

apptainer run  /scratch/anw/share-np/AFNIr MBA -prefix perc_improv \
-chains 6 -iterations 4000 \
-dataTable ${workdir}/MBA_input_acq-${modality}_CORE_upd.txt \
-cVars 'sex' -qVars 'perc_improv,age' \
-stdz 'age,perc_improv' \
-model '1+sex+age+perc_improv' -EOI 'perc_improv' -ridgePlot 40 30  \
-distY 'student' \
-fullRes \
-WCP 4 \
-verb 1

# - responder vs non-responders
modelID=${modality}_responder

mkdir -p ${workdir}/${modelID}
cd ${workdir}/${modelID}

apptainer run  /scratch/anw/share-np/AFNIr MBA -prefix responder \
-chains 6 -iterations 4000 \
-dataTable ${workdir}/MBA_input_acq-${modality}_CORE_upd.txt \
-cVars 'sex,responder' -qVars 'age' \
-stdz 'age' \
-model '1+sex+age+responder' -EOI 'responder' -ridgePlot 40 30  \
-distY 'student' \
-fullRes \
-WCP 4 \
-verb 1


done
