#!/bin/bash
#SBATCH --job-name=MyJob
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --account=williams
#SBATCH --mem=120G
#SBATCH --array=1-40%7

source ~/.bash_profile

homebase=/home/williams/karoraw1
libs_avail=$homebase/2_salmon_inputs/libsheet.txt
lib_prefix=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $libs_avail)
lib_loc=$homebase/transcriptomes

# deal with bzipped read pairs 
lib_type_1_alt=${lib_loc}/${lib_prefix}.fastq.bz2
lib_type_2f_alt=${lib_loc}/${lib_prefix}_1.fastq.bz2
lib_type_2r_alt=${lib_loc}/${lib_prefix}_2.fastq.bz2

if [ -f $lib_type_1_alt ]; then
    bzip2 -d $lib_type_1_alt
    gzip ${lib_loc}/${lib_prefix}.fastq
fi

if [ -f $lib_type_2f_alt ]; then
    bzip2 -d $lib_type_2f_alt
    bzip2 -d $lib_type_2r_alt
    gzip ${lib_loc}/${lib_prefix}_1.fastq
    gzip ${lib_loc}/${lib_prefix}_2.fastq
fi

# return to default gzipped pairs
lib_type_1=${lib_loc}/${lib_prefix}.fastq.gz
lib_type_2f=${lib_loc}/${lib_prefix}_1.fastq.gz
lib_type_2r=${lib_loc}/${lib_prefix}_2.fastq.gz

# get the new index base 
index_base=$homebase/2_salmon_inputs/pulled_genes

if [ -f $lib_type_1 ]; then

    trimLib=${lib_loc}/${lib_prefix}/${lib_prefix}_trimmed.fq.gz    
    if [ ! -f $trimLib ]; then
        conda activate qc
        trim_galore -o $lib_loc/$lib_prefix -j 4 --clip_R1 10 --gzip $lib_type_1
        conda deactivate
    fi

    for index_file in `ls -d $index_base/*idx`; do
        bin_prefix=$(basename $index_file .genes_idx)
        outfile=$homebase/salmon_out2/${lib_prefix}_${bin_prefix}.quant
        salmon quant --mimicBT2 --validateMappings --incompatPrior 0.0 --gcBias -i $index_file --libType IU -r $trimLib -o $outfile --meta -p 40
    done
fi

if [ -f $lib_type_2f ]; then
    fwd=$lib_loc/${lib_prefix}/${lib_prefix}_1_val_1.fq.gz 
    rev=$lib_loc/${lib_prefix}/${lib_prefix}_2_val_2.fq.gz 
    if [ ! -f $fwd ]; then
        conda activate qc
        trim_galore -o $lib_loc/$lib_prefix -j 4 --clip_R1 10 --clip_R2 5 --gzip --paired $lib_type_2f $lib_type_2r
        conda deactivate
    fi

    for index_file in `ls -d $index_base/*idx`; do
        bin_prefix=$(basename $index_file .genes_idx)
        outfile=$homebase/salmon_out2/${lib_prefix}_${bin_prefix}.quant
        salmon quant --mimicBT2 --validateMappings --incompatPrior 0.0 -i $index_file --libType IU -1 $fwd -2 $rev -o $outfile --meta -p 40
    done
fi
