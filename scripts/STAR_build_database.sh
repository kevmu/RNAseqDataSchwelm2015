# /usr/bin/env bash

#BSUB -J STAR_build_database_batch_run
#BSUB -n 28
#BSUB -R "span[hosts=1]"
#BSUB -W 60:00
#BSUB -o STAR_build_database_batch_run%J.out
#BSUB -e STAR_build_database_batch_run%J.err

source ~/.bashrc

conda activate STAR_env

# STAR build genome index file.
STAR --runThreadN 28 --genomeSAindexNbases 11 --runMode genomeGenerate --genomeDir /gpfs/home/kmuirhead/pb/STAR_genome_index --genomeFastaFiles /gpfs/home/kmuirhead/pb/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.fna 





