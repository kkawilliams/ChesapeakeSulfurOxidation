#!/bin/bash

#SBATCH
#SBATCH --job-name=cb_qc
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --partition=gpu
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=../Logs/cb_qc.err
#SBATCH --output=../Logs/cb_qc.out

source activate metawrap2-env
base_dir=/home-3/karoraw1@jhu.edu/scratch/Processed_data_group
D_D=$base_dir/CB33_Summer17_NexteraXT/Raw_Unzipped
T_D=$base_dir/CB33_Summer17_NexteraXT/QCd
sample_names=$base_dir/CB33_Summer17_NexteraXT/CB33_Summer17_NexteraXT.samplenames.tsv

mkdir $T_D

while read Sname F_read R_read; do
    metaWRAP read_qc -t 24 -1 $D_D/$F_read -2 $D_D/$R_read -o $T_D/$Sname
done < $sample_names
