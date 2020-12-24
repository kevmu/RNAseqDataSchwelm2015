# /usr/bin/env bash

#BSUB -J fastqc_batch_run
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 60:00
#BSUB -o fastqc_batch_run%J.out
#BSUB -e fastqc_batch_run%J.err

source ~/.bashrc

conda activate fastqc_env

for i in $(ls | grep "ERR"); 
do 
	echo $i; 
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	fastqc -o /gpfs/home/kmuirhead/pb/$i -f fastq /gpfs/home/kmuirhead/pb/$i/$filename\.1\.fastq\.gz /gpfs/home/kmuirhead/pb/$i/$filename\.2\.fastq\.gz 
done

