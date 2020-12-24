

# Download the Schwelm 2015 RNA-seq dataset.

for i in `cat Pb_eH_SraRunInfo.csv | tail -n +2`; do echo $i; sra_id=$(echo $i | cut -d ',' -f1); echo $sra_id; fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $sra_id; done

