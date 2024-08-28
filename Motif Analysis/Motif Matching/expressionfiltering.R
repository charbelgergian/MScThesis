library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)
library(patchwork)
library(universalmotif)
library(jsonlite)
library(ggplot2)
library(ggvenn)

#Loading HOCOMOCOv12 metadata

hocomoco <- flatten(stream_in(file('/rds/general/user/cg2723/home/tfFinal/H12CORE_annotation.jsonl')))
hocomoco <- hocomoco[,c('name','masterlist_info.species.HUMAN.gene_symbol','tf','masterlist_info.tfclass_family')]

#Function that takes output file of tomtom, list of acetylated promoters and motif database metadata
overlapAceHomer <- function(tomtom,ace,hocomoco, q = 10e-3,sample, metric = 'logCPM'){
  tomtom$Sample <- sample 
  names(tomtom)[names(tomtom) == 'q-value'] <- 'q_value'
  names(tomtom)[names(tomtom) == 'p-value'] <- 'p_value'
  tomtom <- left_join(tomtom,hocomoco,join_by(Target_ID == name),suffix = c("", ""))
  
  top <- tomtom %>% subset(q_value < q)
  top <- top %>% group_by(Target_ID) %>% top_n(-1,q_value)
  
  top <- inner_join(top,ace,join_by("masterlist_info.species.HUMAN.gene_symbol"=='SYMBOL'))
  
  top <- top %>% group_by(tf) %>%  top_n(1,get(metric)) %>% top_n(-1,p_value)
  
}
#Like function above, but separate for ease of adjusting and use with DeepSTARR results
overlapAceDeep <- function(tomtom,ace,hocomoco, q = 10e-3,sample,metric = 'logCPM'){
  tomtom$Sample <- sample
  names(tomtom)[names(tomtom) == 'q-value'] <- 'q_value'
  tomtom <- left_join(tomtom,hocomoco,join_by(Target_ID == name),suffix = c("", ""))
  top <- tomtom %>% subset(q_value < q)
  top <- top %>% group_by(Target_ID) %>% top_n(-1,q_value)
  
  top <- right_join(top,ace,join_by("masterlist_info.species.HUMAN.gene_symbol"=='SYMBOL'))
  top <- top %>% group_by(tf) %>%  top_n(1,get(metric)) %>% top_n(1,q_value)
}

#Function to read output of HOMER and filter for p value
read_filter_homer <- function(homer, p = 1e-12){
  i <- 1
  motifs <- read_homer(homer)
  df <- data.frame()
  for (i in 1:length(motifs)){
    df[i,'name'] <- motifs[[i]]@name
    df[i,'homerPval'] <- motifs[[i]]@pval
  }
  df[df$homerPval<p,]
}

#Reading in output of tomtom for DeepSTARR results
ttDeepERG <- read_delim('/rds/general/user/cg2723/home/tfFinal/tomtom/Deep_ERG/tomtom.tsv',comment = '#')
ttDeepPU1 <- read_delim('/rds/general/user/cg2723/home/tfFinal/tomtom/Deep_PU1/tomtom.tsv',comment = '#')

#Reading in HOMER output and tomtom output
#read_filter_homer keeps only the motifs from HOMER with p-values below the given threshold 
#Inner join is then used to only keep matches in tomtom from those motifs

homerPU1 <- read_filter_homer('/rds/general/user/cg2723/home/tfFinal/Homer/allpu1_auto_bg_250/homerMotifs.all.motifs')
ttHomerPU1 <- read_delim('/rds/general/user/cg2723/home/tfFinal/tomtom/homer_pu1_auto_bg_250/tomtom.tsv',comment = '#')
ttHomerPU1 <- inner_join(ttHomerPU1,homerPU1,join_by("Query_ID" == "name"))

homerERG <- read_filter_homer('/rds/general/user/cg2723/home/tfFinal/Homer/allerg_auto_bg_250/homerMotifs.all.motifs')
ttHomerERG <- read_delim('/rds/general/user/cg2723/home/tfFinal/tomtom/homer_erg_auto_bg_250/tomtom.tsv',comment = '#')
ttHomerERG <- inner_join(ttHomerERG,homerERG,join_by("Query_ID" == "name"))

#Reading in the acetylation data
pu1Diff <- read_delim('/rds/general/user/cg2723/home/tfFinal/Diff/ADvas_AAA_20240422_ac_hg19.pu1VSall.0.01.promoters.txt')
ergDiff <- read_delim('/rds/general/user/cg2723/home/tfFinal/Diff/ADvas_AAA_20240422_ac_hg19.ergVSall.0.01.promoters.txt')

#Overlap tomtom matches with acetylation data
TopDeepERG <- overlapAceDeep(ttDeepERG,ergDiff,hocomoco,sample ='ERG')
TopDeepPU1 <- overlapAceDeep(ttDeepPU1,pu1Diff,hocomoco,sample = 'PU1')

TopHomerERG <- overlapAceHomer(ttHomerERG,ergDiff,hocomoco,sample = 'ERG')
TopHomerPU1 <- overlapAceHomer(ttHomerPU1,pu1Diff,hocomoco,sample = 'PU1')

#Example of heatmap generation 
pheatmap(data.matrix(TopHomerERG[,c('logCPM')],
                     rownames.force = TRUE),
         cluster_rows = F, 
         cluster_cols = F,
         main = '    ERG Homer Auto BG CPM', 
         filename = '/rds/general/user/cg2723/home/tfFinal/plots/ERG_Homer_Auto_250_CPM.png',
         labels_row = TopHomerERG$masterlist_info.species.HUMAN.gene_symbol,
         height = 8, width = 8,
         annotation_row = data.frame(annotationrowHomerERG))
