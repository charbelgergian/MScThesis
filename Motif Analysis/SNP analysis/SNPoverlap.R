library(tidyverse)
library(readxl)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(dplyr)
library(stringr)

hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
ergTF <-import('/rds/general/user/cg2723/home/tfFinal/erg_ADvas_AAA_20240422_tf.bed')
ergAC <-import('/rds/general/user/cg2723/home/tfFinal/erg_ADvas_AAA_20240422_ac.bed')
pu1TF <- import('/rds/general/user/cg2723/home/tfFinal/pu1_ADvas_AAA_20240422_tf.bed')
pu1AC <-import('/rds/general/user/cg2723/home/tfFinal/pu1_ADvas_AAA_20240422_ac.bed')

#Reading in the output from fimo

HomerERGFimo <- read_delim('/rds/general/user/cg2723/home/tfFinal/fimo/HomerERG/fimo.tsv')
HomerERGFimo$motif_id <- str_split_i(HomerERGFimo$motif_id,'\\.',1)

DeepERGFimo <- read_delim('/rds/general/user/cg2723/home/tfFinal/fimo/DeepERG/fimo.tsv')
DeepERGFimo$motif_id <- str_split_i(DeepERGFimo$motif_id,'\\.',1)

HomerPU1Fimo <- read_delim('/rds/general/user/cg2723/home/tfFinal/fimo/HomerPU1/fimo.tsv')
HomerPU1Fimo$motif_id <- str_split_i(HomerPU1Fimo$motif_id,'\\.',1)

DeepPU1Fimo <- read_delim('/rds/general/user/cg2723/home/tfFinal/fimo/DeepPU1/fimo.tsv')
DeepPU1Fimo$motif_id <- str_split_i(DeepPU1Fimo$motif_id,'\\.',1)

#Importing summary stats of GWAS and converting to GRanges
jansenSS <- read_delim('/rds/general/user/cg2723/home/tfFinal/GWAS/AD_Jansen2019_munged.GRCh38.tsv')
jansenSS <- jansenSS[jansenSS$P < 5e-8,]
jansenSS$CHR <- paste('chr',jansenSS$CHR,sep='')
jss <- GRanges(seqnames = jansenSS$CHR, ranges = IRanges(start = jansenSS$BP,width = 1), SNP = jansenSS$SNP, pval = jansenSS$P)

WMHss <- read_tsv('/rds/general/user/cg2723/home/tfFinal/GWAS/SVD_Sargurupremraj2022_munged.GRCh38.tsv')
WMHss <- WMHss[WMHss$P < 5e-8,]
WMHss$CHR <- paste('chr',WMHss$CHR,sep='')
WMHss <- GRanges(seqnames = WMHss$CHR, ranges = IRanges(start = WMHss$BP,width = 1), SNP = WMHss$SNP, pval = WMHss$P)

HomerERGFimo <- drop_na(HomerERGFimo,start,stop)
DeepERGFimo <- drop_na(DeepERGFimo,start,stop)

HomerPU1Fimo <- drop_na(HomerPU1Fimo,start,stop)
DeepPU1Fimo <- drop_na(DeepPU1Fimo,start,stop)

#Converting the fimo results into GRanges objects which contain locations of motifs

HomerERG_GR <- GRanges(seqnames = HomerERGFimo$sequence_name, ranges = IRanges(start = HomerERGFimo$start,end = HomerERGFimo$stop),
                       strand = HomerERGFimo$strand,
                       motif = HomerERGFimo$motif_id,
                       fimostart = HomerERGFimo$start,
                       fimostop = HomerERGFimo$stop)


DeepERG_GR <- GRanges(seqnames = DeepERGFimo$sequence_name, ranges = IRanges(start = DeepERGFimo$start,end = DeepERGFimo$stop),
                      motif = DeepERGFimo$motif_id,
                      fimostart = DeepERGFimo$start,
                      fimostop = DeepERGFimo$stop)

HomerPU1_GR <- GRanges(seqnames = HomerPU1Fimo$sequence_name, ranges = IRanges(start = HomerPU1Fimo$start,end = HomerPU1Fimo$stop),
                       motif = HomerPU1Fimo$motif_id,
                       fimostart = HomerPU1Fimo$start,
                       fimostop = HomerPU1Fimo$stop )

DeepPU1_GR <- GRanges(seqnames = DeepPU1Fimo$sequence_name, ranges = IRanges(start = DeepPU1Fimo$start,end = DeepPU1Fimo$stop),
                      motif = DeepPU1Fimo$motif_id,
                      fimostart = DeepPU1Fimo$start,
                      fimostop = DeepPU1Fimo$stop)

#Overlapping the SNPs with the motif locations
GWAShits <- function(GWAS,gr,ac){
  hits <- join_overlap_intersect(GWAS,gr)
  hits <- subsetByOverlaps(hits,ac)
  hits <- as.GRanges(annotatePeak(hits,TxDb =hg38, annoDb = "org.Hs.eg.db", tssRegion = c(-2000,500)))
}

hitsERGHomer_AD <- GWAShits(jss,HomerERG_GR, ergAC) 
hitsERGDeep_AD <- GWAShits(jss,DeepERG_GR,ergAC)
hitsPU1Homer_AD <- GWAShits(jss,HomerPU1_GR,pu1AC)
hitsPU1Deep_AD <- GWAShits(jss,DeepPU1_GR,pu1AC)

hitsERGHomer_SVD <- GWAShits(WMHss,HomerERG_GR, ergAC)
hitsERGDeep_SVD <- GWAShits(WMHss,DeepERG_GR,ergAC)
hitsPU1Homer_SVD <- GWAShits(WMHss,HomerPU1_GR,pu1AC)
hitsPU1Deep_SVD <- GWAShits(WMHss,DeepPU1_GR,pu1AC)

