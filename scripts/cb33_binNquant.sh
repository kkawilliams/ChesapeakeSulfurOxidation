#!/bin/bash

#SBATCH
#SBATCH --job-name=CB33
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=117600MB
#SBATCH --exclusive
#SBATCH --partition=shared
#SBATCH --mail-type=END
#SBATCH --mail-user=k.arorawilliams2@gmail.com
#SBATCH --error=cb3_bins.err
#SBATCH --output=cb3_bins.out

source activate metawrap2-env

RAW_DIR=/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/Processed_data_group
MW_O=/home-3/karoraw1@jhu.edu/scratch/ChesBayTransect/sandbox
NAME=CB33
READS_DIR=$RAW_DIR/CB33_Summer17_NexteraXT/QCd

# make a single location for relavent libraries
mkdir -p $READS_DIR/QC_Renamed_Seqs

# link them into place
for i in `ls -d $READS_DIR/*`; do
 pfx=`basename $i`
 ln $i/final_pure_reads_1.fastq $READS_DIR/QC_Renamed_Seqs/${pfx}_1.fastq
 ln $i/final_pure_reads_2.fastq $READS_DIR/QC_Renamed_Seqs/${pfx}_2.fastq
done;

# bin them locally 
BIN_HOME=$MW_O/${NAME}_Bins
ASSEMBLY=$MW_O/minimized_contigs/${NAME}_clean_cntgs.fa
N_CPUS=24

metaWRAP binning -t $N_CPUS --maxbin2 -a $ASSEMBLY -o $BIN_HOME $READS_DIR/QC_Renamed_Seqs/*.fastq

# quantify them bins
BINS_DIR=$BIN_HOME/maxbin2_bins

metawrap quant_bins \
-t $N_CPUS \
-b $BINS_DIR \
-o $BIN_HOME/Bin_Quant \
-a $ASSEMBLY \
$READS_DIR/QC_Renamed_Seqs/*fastq

