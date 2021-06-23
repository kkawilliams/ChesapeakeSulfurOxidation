#!/bin/bash

#SBATCH
#SBATCH --job-name=gtdbtk_bins
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --partition=lrgmem
#SBATCH --mem=900G
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --output=gtdbtk_bins.out
#SBATCH --error=gtdbtk_bins.err


bindir=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/data/all_bins
gtdb_tk_out=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/data/gtdb_taxout

ml python/3.7.4-anaconda
source activate gtdbtk
echo `which python`
gtdbtk classify_wf --force -x fa --genome_dir ${bindir} --out_dir $gtdb_tk_out --cpus 48
