#!/bin/bash

#SBATCH
#SBATCH --job-name=@SID@_demux
#SBATCH --time=4:00:00
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --partition=parallel
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=@SID@_demux_test.err
#SBATCH --output=@SID@_demux_test.out

# these will be filled in automatically
SCRIPTS_=~/scratch/ChesBayTransect/scripts
BASE_OUT=/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/Processed_data_group/@SID@
SEQ_ID=FASTQ
IDX_PATH=$BASE_OUT/$SEQ_ID/Undetermined_S0_L001_I1_001.fastq
FWD_PATH=$BASE_OUT/$SEQ_ID/Undetermined_S0_L001_R1_001.fastq
REV_PATH=$BASE_OUT/$SEQ_ID/Undetermined_S0_L001_R2_001.fastq
BCODE=$BASE_OUT/@SID@_barcodes.txt

gzip -d $BASE_OUT/$SEQ_ID/*.gz

source activate base2
python $SCRIPTS_/check_headers.py $PP_DIR
python $SCRIPTS_/rev_comp_bcodes.py $BCODE

# auto path shortcuts
DEMUX_DIR=$BASE_OUT/$SEQ_ID/Demux
#mkdir -p $DEMUX_DIR

echo "Parsing Index File"
ml parallel
echo `date`
parallel -j 24 -a $BCODE.rc --colsep '\t' "grep --no-group-separator -B1 {2} $IDX_PATH | grep "^@" | sed 's/^@//g' > $DEMUX_DIR/{1}.headers"
echo `date`

echo "Demultiplexing" 
ml java
BB_HOME=~/scratch/ChesBayTransect/bin/bbmap

for header in `ls $DEMUX_DIR/*.headers`; do
    path_name=$(dirname $header)
    file_name=$(basename $header)
    sample_name=$(cut -d "." -f 1 <<< "$file_name")
    sample_path=$path_name/$sample_name
    $BB_HOME/filterbyname.sh in=$FWD_PATH in2=$REV_PATH out=$sample_path.R1.fastq out2=$sample_path.R2.fastq names=$header include=t;
done


