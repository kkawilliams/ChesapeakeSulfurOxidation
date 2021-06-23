#!/bin/bash

#SBATCH
#SBATCH --job-name=dadaEnd
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=960G
#SBATCH --exclusive
#SBATCH --partition=lrgmem
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=../Logs/merge_taxa.err
#SBATCH --output=../Logs/merge_taxa.out

module load R/3.5.1

# these will be filled in automatically
BASE_DIR=/home-3/karoraw1@jhu.edu/scratch/ChesBayTransect
SCRIPTS_=$BASE_DIR/scripts
TAX_DB=$BASE_DIR/data/Silva_DB
LIBdirDF=$BASE_DIR/data/TrimOTUsData/seq_tabs.csv
OUT_DIR=$BASE_DIR/data/TrimOTUsData

Rscript $SCRIPTS_/merge_chim_tax.R $TAX_DB $LIBdirDF $OUT_DIR
