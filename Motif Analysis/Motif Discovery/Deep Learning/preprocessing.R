#!/usr/bin/env Rscript

#This script takes a bed file and a bam file and processes them for the model
#A dataframe is created with three columns: the chromosome number, the DNA sequence and the TF binding score

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)


summits <- import(summitbed)
start(summits) <- start(summits) - 125
end(summits) <- end(summits) + 125

bam <- import(bam)
overlaps <- countOverlaps(summits, bam)

seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38,summits)
seq <- as.list(as.character(seq))
chr <- as.vector(seqnames(summits))
df <- cbind(chr, unlist(seq), overlaps)
df <- data.frame(df)
colnames(df) <- c("chr", "seq", "score")


write.csv(df, outputPath, row.names = FALSE)
