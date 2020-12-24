# /usr/bin/env bash

#BSUB -J trimmomatic_batch_run
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 60:00
#BSUB -o trimmomatic_batch_run%J.out
#BSUB -e trimmomatic_batch_run%J.err

source ~/.bashrc

conda activate trimmoamatic_env

for i in $(ls | grep "ERR"); 
do 
	echo $i;
	mkdir -p "/gpfs/home/kmuirhead/pb/$i/trimmomatic_dir" 
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	trimmomatic PE /gpfs/home/kmuirhead/pb/$i/${filename}.1.fastq.gz /gpfs/home/kmuirhead/pb/$i/${filename}.2.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.paired.1.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.unpaired.1.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.paired.2.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.unpaired.2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

done

