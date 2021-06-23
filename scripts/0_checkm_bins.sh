#!/bin/bash

#SBATCH
#SBATCH --job-name=checkm_bins
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=parallel,shared
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --output=checkm_bins.out
#SBATCH --error=checkm_bins.err


ml python/3.7.4-anaconda
source activate checkenv

bindir=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/data/all_bins

checkmout=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/data/checkmout

tmpdir=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/data/checkm_tmp

export PATH=$PATH:~/.local/lib/python3.7/site-packages/checkm:~/.local/bin

mkdir $tmpdir
checkm lineage_wf --tmpdir $tmpdir -t 24 -x fa $bindir $checkmout
rm -rf $tmpdir
