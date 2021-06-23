
datadir=/Volumes/KeithSSD/SulFox/data/possible_data_sources
tableFile=$datadir/SraRunTable.txt

fastqdir=$datadir/$(cut -d, -f6 $tableFile | head -2 | tail -1)
mkdir -p $fastqdir

for SRR_ID in $(cut -d, -f1 $tableFile | tail -5); do 
    fastq-dump -v --outdir $fastqdir --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $SRR_ID
done


