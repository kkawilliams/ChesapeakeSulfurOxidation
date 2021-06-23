#!/bin/bash -l

#SBATCH

#SBATCH --job-name=*SN*_asm
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100GB
#SBATCH --exclusive
#SBATCH --partition=parallel
#SBATCH --mail-type=END
#SBATCH --mail-user=k.arorawilliams2@gmail.com
#SBATCH --error=*SN*_asm.err
#SBATCH --output=*SN*_asm.out

source activate base2

PE1=*FWD* # full path
PE2=*REV* # full path
OUT=*OUT* # full path
SAMP=*SN* # name

mkdir -p $OUT
mkdir -p $OUT/${SAMP}_tmp

megahit -t 24 -1 $PE1 -2 $PE2 -o $OUT/$SAMP --min-contig-len 500 --tmp-dir $OUT/${SAMP}_tmp -m 0.5 --continue

rm -rf $OUT/${SAMP}_tmp
