#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(cummeRbund)

cuff <- readCufflinks("/mnt/d/cuffdiff")

myGeneIds <- c("PBRA_005464","PBRA_005765","PBRA_003743","PBRA_008059","PBRA_007378","PBRA_003759","PBRA_001907","PBRA_002543","PBRA_002958","PBRA_007091","PBRA_002230","PBRA_008962","PBRA_006006","PBRA_003161","PBRA_005081","PBRA_002551","PBRA_007776","PBRA_001295","PBRA_008942","PBRA_004239")

myGenes <- getGenes(cuff, myGeneIds)

pdf(file="chitin_heatmap.pdf")


csHeatmap(myGenes,cluster='both')
dev.off()
