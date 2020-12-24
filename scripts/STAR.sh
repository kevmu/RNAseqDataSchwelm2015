# /usr/bin/env bash

#BSUB -J STAR_batch_run
#BSUB -n 56
#BSUB -W 60:00
#BSUB -R "rusage[mem=256000]"
#BSUB -o STAR_batch_run%J.out
#BSUB -e STAR_batch_run%J.err

source ~/.bashrc

conda activate STAR_env

for i in $(ls | grep "ERR"); 
do 
	echo $i;
	STAR_dir="/gpfs/home/kmuirhead/pb/$i/STAR_dir" 
	mkdir -p $STAR_dir
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	
	gzipped_fastq_file1="/gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.paired.1.fastq.gz"
	gzipped_fastq_file2="/gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.paired.2.fastq.gz"

	fastq_file1=$(echo $gzipped_fastq_file1 | sed 's/\.gz//g')
	fastq_file2=$(echo $gzipped_fastq_file2 | sed 's/\.gz//g')

	echo $fastq_file1
	echo $fastq_file2

	gunzip -c $gzipped_fastq_file1 > $fastq_file1
	gunzip -c $gzipped_fastq_file2 > $fastq_file2

	#exit 0;
	# STAR mapping for scRNA-seq pbrassicae.
	STAR --runThreadN 56 --runMode alignReads --genomeDir /gpfs/home/kmuirhead/pb/STAR_genome_index --readFilesIn $fastq_file1 $fastq_file2 --sjdbGTFfile /gpfs/home/kmuirhead/pb/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.gtf --outFileNamePrefix $STAR_dir/${filename}_ --outSAMtype SAM --outSAMattributes All --outSAMstrandField intronMotif


done

