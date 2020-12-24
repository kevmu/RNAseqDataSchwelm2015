Kevin Muirhead and Edel Pérez-López, "Plasmodiophora brassicae chitin-binding effectors guard and mask spores during infection" 

## RNA-Seq Data Analysis Workflow (Schwelle 2015)

# Download the Schwelm 2015 RNA-seq dataset.
for i in `cat Pb_eH_SraRunInfo.csv | tail -n +2`; do echo $i; sra_id=$(echo $i | cut -d ',' -f1); echo $sra_id; fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $sra_id; done

# Use md5sum command to check the integrity of the md5 checksum files downloaded with the SRA dataset.
for i in $(ls | grep "ERR"); do echo $i; cd $i; md5sum -c $i.md5 >> md5checksum_results.txt; cd /gpfs/home/kmuirhead/pb; done
