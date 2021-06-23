#!/bin/bash

#SBATCH
#SBATCH --job-name=decompress
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=parallel,shared
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --output=decompress.out
#SBATCH --error=decompress.err

mappingFile=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/data/readmapping/batches/BATCH0.txt

ml samtools parallel
mapHELP=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/scripts/mapHELPER.sh

parallel --colsep '\t' "mapHELP {1} {2} {3} {4}" :::: $mappingFile
