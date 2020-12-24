# /usr/bin/env bash

#BSUB -J samtools_batch_run
#BSUB -n 56
#BSUB -W 60:00
#BSUB -R "rusage[mem=256000]"
#BSUB -o samtools_batch_run%J.out
#BSUB -e samtools_batch_run%J.err

source ~/.bashrc

conda activate samtools_env

#for i in $(ls | grep "ERR");
for i in $(ls | grep "ERR1337811");  
do 
	echo $i;
	STAR_dir="/gpfs/home/kmuirhead/pb/$i/STAR_dir"
	sam_file=$(ls $STAR_dir/* | grep "sam");
	sam_filename=$(basename $sam_file);
	bam_filename=$(echo $sam_filename | sed 's/\.sam/\.bam/g'); 
	samtools sort -O bam -o $STAR_dir/$bam_filename -@ 56 $STAR_dir/$sam_filename; 
done

