r1=$1
r2=$2
combo=$3
prefix=$4

if [ ! -f ${combo}.depth ]; then

    BWA2_exe=/home-3/karoraw1@jhu.edu/scratch/CBFunctions/bin/bwa-mem2-2.0pre2_x64-linux/bwa-mem2
    $BWA2_exe mem -t 1 ${prefix} ${r1} ${r2} > ${combo}.sam

    # convert
    ml samtools
    samtools view -@ 1 -bS ${combo}.sam > ${combo}.bam

    # depth
    samtools depth -d 300 -a ${combo}.bam > ${combo}.depth

fi 
# delete
rm ${combo}.sam
rm ${combo}.bam
