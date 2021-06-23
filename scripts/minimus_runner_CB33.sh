#!/bin/bash

#SBATCH
#SBATCH --job-name=miniCB3
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --partition=lrgmem
#SBATCH --mem=800G
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=../Logs/miniCB3.err
#SBATCH --output=../Logs/miniCB3.out

source activate amosEnv
PERL5LIB=/scratch/users/karoraw1@jhu.edu/miniconda2/envs/amosEnv/lib
PERL5LIB=${PERL5LIB}:/scratch/users/karoraw1@jhu.edu/miniconda2/envs/amosEnv/lib/TIGR
PERL5LIB=${PERL5LIB}:/scratch/users/karoraw1@jhu.edu/miniconda2/envs/amosEnv/lib/AMOS
export PERL5LIB=$PERL5LIB

home_dir=$(dirname $(pwd))
assem_list=${home_dir}/data/cb33_assemblies.txt
agg_dir=`pwd`/Minimus_Out_CB33
mini_conts=$agg_dir/CB33_miniConts_plusSings.fa
mini_conts2=$agg_dir/CB33_dd_contigs.fa

mkdir -p $agg_dir

temp_file=$agg_dir/edited_contigs.temp.fasta
agg_file=$agg_dir/contig_pile.fasta
afg_file=$agg_dir/minimized_cntgs.afg
cntg_file=${afg_file%.*}.fasta
singletons=${afg_file%.*}.singletons.seq

while read assem_file; do
 pfx=$(basename $(dirname ${assem_file}))  
 sed "s/>k141/>${pfx}_k141/g" $assem_file > ${temp_file}
 cat ${temp_file} >> ${agg_file}
 echo "Contigs in agg file = " $(grep -c ">" ${agg_file}) 
done < $assem_list

rm $temp_file
ml java
$home_dir/bin/bbmap/dedupe.sh in=${agg_file} out=${mini_conts2} maxedits=50 minidentity=97

toAmos -s $mini_conts2 -o $afg_file
cd $agg_dir
minimus2 $(basename $afg_file .afg) -D OVERLAP=80
cat $cntg_file $singletons > $mini_conts

