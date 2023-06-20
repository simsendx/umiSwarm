#!/bin/bash

# Job name: dada2 bacteria soils
#SBATCH --job-name=dada2__a_full
#SBATCH --account=nn9744k
#SBATCH --time=3:00:00
#SBATCH --mem=61GB
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
##SBATCH --array 1-26

##Setup job environment
module purge   # clear any inherited modules
set -o errexit # exit on errors

module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0


Rscript /cluster/projects/nn9745k/scripts/40_simsenseq/dada2_fish20230313_simsen.R