# /usr/bin/env bash

#BSUB -J cuffdiff_batch_run
#BSUB -n 56
#BSUB -W 60:00
#BSUB -R "rusage[mem=256000]"
#BSUB -o cuffdiff_batch_run%J.out
#BSUB -e cuffdiff_batch_run%J.err

source ~/.bashrc

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


