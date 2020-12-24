Kevin Muirhead and Edel Pérez-López, "Plasmodiophora brassicae chitin-binding effectors guard and mask spores during infection" 


### Supplementary Data 


## NCBI SRA Run Information File. Contains SRA accession IDs and other metadata for the RNA-Seq experiment.

Pb_eH_SraRunInfo.csv



## Scripts

scripts/

fastqc.sh
trimmomatic.sh
STAR_build_database.sh
STAR.sh
samtools.sh
cuffdiff.sh
cummeRbund_heatmap.R


## Conda Environments

envs/

# The NCBI sra-toolkit conda environment.
sra_toolkit_env.yml

# The fastqc condo environment.
fastqc_env.yml

# The trimmoamatic condo environment.
trimmoamatic_env.yml

# The STAR condo environment.
STAR_env.yml

# The samtools condo environment.
samtools_env.yml

# The cufflinks condo environment.
cufflinks_env.yml


## Dropbox Data

# cuffdiff output directory containing files for generating the cummerRbund heatmap.
https://www.dropbox.com/sh/mhjim865vfoohu9/AADGqJ7_wluX5eegvJtb6T7ia?dl=0

# cuffdiff Legend
ID, SRA Accession ID, RNA-Seq Library Name, Heatmap Name, Heatmap Category
q1, ERR1337805, e3Brapa, Br, Plant Host
q2, ERR1337806, e3Bnapus, Bn, Plant Host
q3, ERR1337807, e3Bole, Bo, Plant Host
q4, ERR1337808, e3Germ, Gs, Development
q5, ERR1337809, e3Mature, Ms, Development
q6, ERR1337810, e3Mix, N/A, N/A
q7, ERR1337811, e3Plasmodia, Pl, Development



## RNA-Seq Data Analysis Workflow using the Schwelle 2015 RNA-Seq Dataset

## Download the P. brassicae Pbe3.h15 genome assembly (Schwelm et al. 2015) for mapping RNA-seq reads using STAR (Dobin et al. 2013)

# Link to the Pbe3.h15 genome at the NCBI assembly database.
https://www.ncbi.nlm.nih.gov/assembly/GCA_001049375.1

# Link to the Genbank GCA assembly on the NCBI FTP site.
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/049/375/GCA_001049375.1_pbe3.h15/

# Download Genbank GCA assembly on the NCBI FTP site using the wget command.
wget --recursive --no-host-directories --cut-dirs=6 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/049/375/GCA_001049375.1_pbe3.h15/ -P my_dir/


## If you need to install the conda environments we used in this article. Use the following conda commands for creating the following environments.

# Create the NCBI sra-toolkit condo environment.
conda create env --file sra_toolkit_env.yml

# Create the fastqc condo environment.
conda create env --file fastqc_env.yml

# Create the trimmoamatic condo environment.
conda create env --file trimmoamatic_env.yml

# Create the STAR condo environment.
conda create env --file STAR_env.yml

# Create the samtools condo environment.
conda create env --file samtools_env.yml

# Create the cufflinks condo environment.
conda create env --file cufflinks_env.yml


## Start of RNA-Seq Data Analysis Workflow 

# Create working directory.
working_dir="$working_dir"
mkdir -p $working_dir

# Change directory to working directory.
cd $working_dir

# Activate the NCBI sra-toolkit conda environment.
conda activate sra_toolkit_env

# Download the Schwelm 2015 RNA-seq SRA dataset using the fastq-dump program of the the NCBI sra-toolkit.
for i in $(cat "${working_dir}/Pb_eH_SraRunInfo.csv" | tail -n +2); do echo $i; sra_id=$(echo $i | cut -d ',' -f1); echo $sra_id; fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $sra_id; done

# Change directory to working directory.
cd $working_dir

# Use md5sum command to check the integrity of the md5 checksum files downloaded with the SRA dataset.
for i in $(ls | grep "ERR"); do echo $i; cd $i; md5sum -c $i.md5 >> $working_dir/md5checksum_results.txt; cd $working_dir; done

# Rename the fastq files in the ERR1337805 to the same suffix as the other RNA-Seq libraries so that we can use a for loop to iterate through each library.
cd $working_dir
cd ERR1337805
cp E3BR_l1_1.fq.gz E3BR_l1.1.fastq.gz
cp E3BR_l1_2.fq.gz E3BR_l1.2.fastq.gz


### FastQC for quality control assessment and to figure out which sequencing adapters were used.

# Activate the FastQC conda environment.
conda activate fastqc_env

# Iterate through each RNA-seq Library directory and perform quality control assessment using the fastqc program.
for i in $(ls | grep "ERR"); 
do 
	echo $i; 
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	fastqc -o $working_dir/$i -f fastq $working_dir/$i/$filename\.1\.fastq\.gz $working_dir/$i/$filename\.2\.fastq\.gz 
done


# Activate the trimmoamatic conda environment.
conda activate trimmoamatic_env

# Iterate through each RNA-seq Library directory and perform quality-filtering and trimming of adapters using trimmomatic.
for i in $(ls | grep "ERR"); 
do 
	echo $i;
	mkdir -p "$working_dir/$i/trimmomatic_dir" 
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	trimmomatic PE $working_dir/$i/${filename}.1.fastq.gz $working_dir/$i/${filename}.2.fastq.gz $working_dir/$i/trimmomatic_dir/${filename}.paired.1.fastq.gz $working_dir/$i/trimmomatic_dir/${filename}.unpaired.1.fastq.gz $working_dir/$i/trimmomatic_dir/${filename}.paired.2.fastq.gz $working_dir/$i/trimmomatic_dir/${filename}.unpaired.2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

done

### STAR Mapping

# Activate the STAR conda environment.
conda activate STAR_env

# STAR build genome index file using 28 cpu threads. (scripts/STAR_build_database.sh)
STAR --runThreadN 28 --genomeSAindexNbases 11 --runMode genomeGenerate --genomeDir $working_dir/STAR_genome_index --genomeFastaFiles $working_dir/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.fna 


# Mapping RNA-seq reads to the P. brassicae Pbe3.h15 genome (Schwelm et al. 2015) using STAR. Using 56 cpu threads. (scripts/STAR.sh)
for i in $(ls | grep "ERR"); 
do 
	echo $i;
	STAR_dir="$working_dir/$i/STAR_dir" 
	mkdir -p $STAR_dir
	filename=$(ls $i | grep ".gz" | sed 's/\.[1-2]\.fastq\.gz//g' | sort | uniq); 
	echo $filename;
	
	gzipped_fastq_file1="$working_dir/$i/trimmomatic_dir/${filename}.paired.1.fastq.gz"
	gzipped_fastq_file2="$working_dir/$i/trimmomatic_dir/${filename}.paired.2.fastq.gz"

	fastq_file1=$(echo $gzipped_fastq_file1 | sed 's/\.gz//g')
	fastq_file2=$(echo $gzipped_fastq_file2 | sed 's/\.gz//g')

	echo $fastq_file1
	echo $fastq_file2

	gunzip -c $gzipped_fastq_file1 > $fastq_file1
	gunzip -c $gzipped_fastq_file2 > $fastq_file2

	# STAR mapping for scRNA-seq pbrassicae.
	STAR --runThreadN 56 --runMode alignReads --genomeDir $working_dir/STAR_genome_index --readFilesIn $fastq_file1 $fastq_file2 --sjdbGTFfile $working_dir/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.gtf --outFileNamePrefix $STAR_dir/${filename}_ --outSAMtype SAM --outSAMattributes All --outSAMstrandField intronMotif


done


### Sorting SAM to BAM using samtools.

# Activate the samtools conda environment.
conda activate samtools_env

# Samtools for sorting SAM format and converting to BAM format using 56 cpu threads (scripts/samtools.sh).
for i in $(ls | grep "ERR");
do 
	echo $i;
	STAR_dir="$working_dir/$i/STAR_dir"
	sam_file=$(ls $STAR_dir/* | grep "sam");
	sam_filename=$(basename $sam_file);
	bam_filename=$(echo $sam_filename | sed 's/\.sam/\.bam/g'); 
	samtools sort -O bam -o $STAR_dir/$bam_filename -@ 56 $STAR_dir/$sam_filename; 
done



### Differential expression analysis using cuffdiff program from cufflinks.


# Activate the cufflinks conda environment.
conda activate cufflinks_env

# Perform differential expression analysis using cuffdiff program from cufflinks (scripts/cufflinks.sh).

cuffdiff $working_dir/STAR_genome_index/GCA_001049375.1_pbe3.h15_genomic.gtf \
$working_dir/ERR1337805/STAR_dir/E3BR_l1_Aligned.out.bam \
$working_dir/ERR1337806/STAR_dir/116_dual75_Aligned.out.bam \
$working_dir/ERR1337807/STAR_dir/117_dual76_Aligned.out.bam \
$working_dir/ERR1337808/STAR_dir/105F_index25_Aligned.out.bam \
$working_dir/ERR1337809/STAR_dir/103F_index22_Aligned.out.bam \
$working_dir/ERR1337810/STAR_dir/102F_index21_Aligned.out.bam \
$working_dir/ERR1337811/STAR_dir/P820_101_Aligned.out.bam \
-p 56 -o $working_dir/cuffdiff


# CummeRbund for generating a heatmap to visualize the differential expression analysis from cuffdiff. This is an Rscript (scripts/cummeRbund_heatmap.R).


#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(cummeRbund)

cuff <- readCufflinks("/mnt/d/cuffdiff")

myGeneIds <- c("PBRA_005464","PBRA_005765","PBRA_003743","PBRA_008059","PBRA_007378","PBRA_003759","PBRA_001907","PBRA_002543","PBRA_002958","PBRA_007091","PBRA_002230","PBRA_008962","PBRA_006006","PBRA_003161","PBRA_005081","PBRA_002551","PBRA_007776","PBRA_001295","PBRA_008942","PBRA_004239")

myGenes <- getGenes(cuff, myGeneIds)

pdf(file="chitin_heatmap.pdf")


csHeatmap(myGenes,cluster='both')
dev.off()





