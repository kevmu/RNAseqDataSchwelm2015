## Kevin Muirhead and Edel Pérez-López, "Plasmodiophora brassicae chitin-binding effectors guard and mask spores during infection" 

## RNA-Seq Data Analysis Workflow (Schwelle 2015)

# Create working directory.
working_dir="/gpfs/home/kmuirhead/pb"
mkdir -p $working_dir

# Change directory to working directory.
cd $working_dir

# Download the Schwelm 2015 RNA-seq dataset.
for i in $(cat "${working_dir}/Pb_eH_SraRunInfo.csv" | tail -n +2); do echo $i; sra_id=$(echo $i | cut -d ',' -f1); echo $sra_id; fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $sra_id; done

# Change directory to working directory.
cd $working_dir

# Use md5sum command to check the integrity of the md5 checksum files downloaded with the SRA dataset.
for i in $(ls | grep "ERR"); do echo $i; cd $i; md5sum -c $i.md5 >> $working_dir/md5checksum_results.txt; cd $working_dir; done

# Rename the fastq files in the ERR1337805 to the same suffix as the other RNA-Seq libraries.
cd $working_dir
cd ERR1337805
cp E3BR_l1_1.fq.gz E3BR_l1.1.fastq.gz
cp E3BR_l1_2.fq.gz E3BR_l1.2.fastq.gz


# Activate the FastQC condo environment.
conda activate fastqc_env

# 
for i in $(ls | grep "ERR"); 
do 
	echo $i; 
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	fastqc -o $working_dir/$i -f fastq $working_dir/$i/$filename\.1\.fastq\.gz $working_dir/$i/$filename\.2\.fastq\.gz 
done


#
conda activate trimmoamatic_env

for i in $(ls | grep "ERR"); 
do 
	echo $i;
	mkdir -p "/gpfs/home/kmuirhead/pb/$i/trimmomatic_dir" 
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	trimmomatic PE /gpfs/home/kmuirhead/pb/$i/${filename}.1.fastq.gz /gpfs/home/kmuirhead/pb/$i/${filename}.2.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.paired.1.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.unpaired.1.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.paired.2.fastq.gz /gpfs/home/kmuirhead/pb/$i/trimmomatic_dir/${filename}.unpaired.2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

done


#
conda activate STAR_env

# STAR build genome index file.
STAR --runThreadN 28 --genomeSAindexNbases 11 --runMode genomeGenerate --genomeDir /gpfs/home/kmuirhead/pb/STAR_genome_index --genomeFastaFiles /gpfs/home/kmuirhead/pb/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.fna 


#

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

	# STAR mapping for scRNA-seq pbrassicae.
	STAR --runThreadN 56 --runMode alignReads --genomeDir /gpfs/home/kmuirhead/pb/STAR_genome_index --readFilesIn $fastq_file1 $fastq_file2 --sjdbGTFfile /gpfs/home/kmuirhead/pb/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.gtf --outFileNamePrefix $STAR_dir/${filename}_ --outSAMtype SAM --outSAMattributes All --outSAMstrandField intronMotif


done


###
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


###

conda activate cufflinks_env

cuffdiff /gpfs/home/kmuirhead/pb/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.gtf \
/gpfs/home/kmuirhead/pb/ERR1337805/STAR_dir/E3BR_l1_Aligned.out.bam \
/gpfs/home/kmuirhead/pb/ERR1337806/STAR_dir/116_dual75_Aligned.out.bam \
/gpfs/home/kmuirhead/pb/ERR1337807/STAR_dir/117_dual76_Aligned.out.bam \
/gpfs/home/kmuirhead/pb/ERR1337808/STAR_dir/105F_index25_Aligned.out.bam \
/gpfs/home/kmuirhead/pb/ERR1337809/STAR_dir/103F_index22_Aligned.out.bam \
/gpfs/home/kmuirhead/pb/ERR1337810/STAR_dir/102F_index21_Aligned.out.bam \
/gpfs/home/kmuirhead/pb/ERR1337811/STAR_dir/P820_101_Aligned.out.bam \
-p 56 -o /gpfs/home/kmuirhead/pb/cuffdiff


## Rscript


#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(cummeRbund)

cuff <- readCufflinks("/mnt/d/cuffdiff")

myGeneIds <- c("PBRA_005464","PBRA_005765","PBRA_003743","PBRA_008059","PBRA_007378","PBRA_003759","PBRA_001907","PBRA_002543","PBRA_002958","PBRA_007091","PBRA_002230","PBRA_008962","PBRA_006006","PBRA_003161","PBRA_005081","PBRA_002551","PBRA_007776","PBRA_001295","PBRA_008942","PBRA_004239")

myGenes <- getGenes(cuff, myGeneIds)

pdf(file="chitin_heatmap.pdf")


csHeatmap(myGenes,cluster='both')
dev.off()





