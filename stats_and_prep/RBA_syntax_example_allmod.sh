#!/bin/bash
## C. Vriend - Amsterdam UMC - Aug '24

# slurm settings
#SBATCH --job-name=RBA
#SBATCH --mem=4G
#SBATCH --partition=luna-cpu-long
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=24
#SBATCH --time=01-0:00:00
#SBATCH --nice=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --output=RBA_%A.out


# running containerized version of Region-Based Analysis Program through Bayesian Multilevel Modeling (RBA) available in AFNI
#https://afni.nimh.nih.gov/pub/dist/doc/program_help/RBA.html
# Chen et al. 2019 paper - https://doi.org/10.1016/j.neuroimage.2019.116320


workdir=/data/anw/anw-work/NP/projects/data_chris/CORE/stats/RBA
export APPTAINER_BIND="/data/anw/anw-work,/scratch"

#############
# multilayer 
#############

modality=multi

for graphmeasure in eigenvector eccentricity BC; do 

# When the values for response
#         are too small or large, it may create a convergence problem for MCMC. To
         # avoid the problem, set a scaling factor so that the range of value is
         # around 1-10. The results will be adjusted back to the orignal scale.

# scale 100 for eigenvector / scale 10 for eccentricity

if [[ ${graphmeasure} == "eigenvector" ]]; then

scaleF=100
elif [[ ${graphmeasure} == "eccentricity" ]]; then
scaleF=10

elif [[ ${graphmeasure} == "BC" ]]; then 

scaleF=10

fi 


# - percentage improvement

modelID=${modality}_perc_improv_${graphmeasure}

mkdir -p ${workdir}/${modelID}
cd ${workdir}/${modelID}

eval "apptainer run  /scratch/anw/share-np/AFNIr RBA -prefix perc_improv \
-chains 6 -iterations 4000 \
-dataTable ${workdir}/RBA_input_acq-multi_CORE_upd.txt \
-cVars 'sex' -qVars 'perc_improv,age' \
-stdz 'age,perc_improv' \
-model '1+sex+age+perc_improv' -EOI 'perc_improv' -ridgePlot 40 30  \
-distY 'student' \
-ROI 'network' \
-Y ${graphmeasure} \
-scale ${scaleF} \
-WCP 4 \
-PDP 2 4 \
-verb 1"



# - responder vs non-responders
outcome=responder
modelID=${modality}_${outcome}_${graphmeasure}
mkdir -p ${workdir}/${modelID}
cd ${workdir}/${modelID}


eval "apptainer run  /scratch/anw/share-np/AFNIr RBA -prefix ${outcome} \
-chains 6 -iterations 4000 \
-dataTable ${workdir}/RBA_input_acq-multi_CORE_upd.txt \
-cVars 'sex,${outcome}' -qVars 'age' \
-stdz 'age' \
-model '1+sex+age+${outcome}' -EOI '${outcome}' -ridgePlot 40 30  \
-distY 'student' \
-ROI 'network' \
-Y ${graphmeasure} \
-scale ${scaleF} \
-WCP 4 \
-PDP 2 4 \
-verb 1"

done


#############
# single layer 
#############

for modality in dwi func; do 


# -  percentage improvement
modelID=${modality}_perc_improv_PC

mkdir -p ${workdir}/${modelID}
cd ${workdir}/${modelID}

eval "apptainer run  /scratch/anw/share-np/AFNIr RBA -prefix perc_improv \
-chains 6 -iterations 4000 \
-dataTable ${workdir}/RBA_PC_acq-${modality}_CORE_upd.txt \
-cVars 'sex' -qVars 'perc_improv,age' \
-stdz 'age,perc_improv' \
-model '1+sex+age+perc_improv' -EOI 'perc_improv' -ridgePlot 40 30  \
-distY 'student' \
-ROI 'network' \
-scale 10 \
-WCP 4 \
-PDP 2 4 \
-verb 1"

# - responder vs non-responders
outcome=responder
modelID=${modality}_responder_PC

mkdir -p ${workdir}/${modelID}
cd ${workdir}/${modelID}


eval "apptainer run  /scratch/anw/share-np/AFNIr RBA -prefix ${outcome} \
-chains 6 -iterations 4000 \
-dataTable ${workdir}/RBA_PC_acq-${modality}_CORE_upd.txt \
-cVars 'sex,${outcome}' -qVars 'age' \
-stdz 'age' \
-model '1+sex+age+${outcome}' -EOI '${outcome}' -ridgePlot 40 30  \
-distY 'student' \
-ROI 'network' \
-scale 10 \
-WCP 4 \
-PDP 2 4 \
-verb 1"


done



