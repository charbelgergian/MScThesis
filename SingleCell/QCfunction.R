library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(dplyr)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

createSignac <- function(path,sampleID,samplename){
  metadata <- read.csv(paste(path,'/singlecell.csv',sep=''),header = TRUE,row.names = 1)
  counts <- Read10X_h5(paste(path,'/filtered_peak_bc_matrix.h5',sep=''))
  frag <-  CreateFragmentObject(paste(path,'/fragments.tsv.gz',sep=''))
  assay <- CreateChromatinAssay(counts = counts,sep = c(":", "-"),fragments = frag)
  
  sampleID <- CreateSeuratObject(counts = assay,assay = "peaks_",project=samplename,meta.data = metadata)
  Annotation(sampleID) <- annotations 
  sampleID <- TSSEnrichment(sampleID, fast = FALSE)
  sampleID$frip <- sampleID$peak_region_fragments / sampleID$passed_filters * 100
  #sampleID <- FRiP(sampleID,assay='peaks_',total.fragments='fragments',col.name = 'frip')
  return(sampleID)
}


path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF125165/outs'
IGF125165 <- createSignac(path,'IGF125165','IGF125165')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF131914/outs'
IGF131914 <- createSignac(path,'IGF131914','IGF131914')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF128020/outs'
IGF128020 <- createSignac(path,'IGF128020','IGF128020')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF133738/outs'
IGF133738 <- createSignac(path,'IGF133738','IGF133738')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF133740/outs'
IGF133740 <- createSignac(path,'IGF133740','IGF133740')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF133741/outs'
IGF133741 <- createSignac(path,'IGF133741','IGF133741')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF133742/outs'
IGF133742 <- createSignac(path,'IGF133742','IGF133742')
path <- '/rds/general/user/cg2723/home/SC/submitfiles/IGF133743/outs'
IGF133743 <- createSignac(path,'IGF133743','IGF133743')

samples <- c(IGF125165,IGF128020,IGF131914,IGF133738,IGF133738,IGF133740,IGF133741,IGF133742,IGF133743)
sampleNames <- c('IGF125165','IGF128020','IGF131914','IGF133738','IGF133738','IGF133740','IGF133741','IGF133742','IGF133743')

QCdf <- data.frame()
for (sample in samples){
  sampleDF <- data.frame(frip = sample@meta.data$frip, fragments = sample@meta.data$passed_filters,TSSE = sample@meta.data$TSS.enrichment, sample = sample@project.name)
  QCdf <- rbind(QCdf,sampleDF)
}

ggplot(QCdf,aes(sample,frip))+geom_violin()+
  geom_hline(yintercept =15,color = 'orange')+
  geom_hline(yintercept =25,color = 'green')+
  scale_x_discrete(labels = 1:8)+
  xlab('Sample')+
  ylab('FriP')+
  ggtitle('c')
ggsave('/rds/general/user/cg2723/home/SC/plotsnew/frip.png', width = 6, height = 3, dpi=1200)

ggplot(QCdf,aes(sample,fragments))+scale_y_log10()+geom_boxplot()+
  scale_x_discrete(labels = 1:8)+
  geom_hline(yintercept = 20000, color = 'blue')+
  ggtitle('b')+
  xlab('Sample')+
  ylab('Unique reads')
ggsave('/rds/general/user/cg2723/home/SC/plotsnew/reads.png', width = 6, height = 3, dpi=1200)


ggplot(QCdf,aes(sample,TSSE))+
  geom_violin()+
  geom_hline(yintercept = 3, color = 'orange') +
  geom_hline(yintercept = 5, color ='green')+
  ggtitle('d')+
  xlab('Sample')+
  ylab('TSS Enrichment')
ggsave('/rds/general/user/cg2723/home/SC/plotsnew/tsse.png', width = 6, height = 3, dpi=1200)

counts <- count(QCdf,sample)
ggplot(counts,aes(x=sample, y=n,fill=n))+geom_bar(stat='identity')+
  scale_x_discrete(labels = 1:8)+
  scale_fill_gradient(low = "pink", high = "blue")+
  xlab('Sample')+
  ylab('Count')+
  ggtitle('a')

ggsave('/rds/general/user/cg2723/home/SC/plotsnew/count.png', width = 6, height = 3, dpi=1200)

ggplot(QCdf,aes(fragments,frip))+geom_point(size=0.1)+facet_wrap(~sample)


Norm <- function(Seu){
  Seu <- RunTFIDF(Seu)
  Seu <- FindTopFeatures(Seu, min.cutoff = 'q0')
  Seu <- RunSVD(Seu)
  return(Seu)
}    

UMAP <- function(Seu,dims=2:30,resolution=1){
  Seu <- RunUMAP(Seu,reduction = 'lsi', dims = dims)
  Seu <- FindNeighbors(Seu,reduction = 'lsi', dims = dims)
  Seu <- FindClusters(Seu, verbose = FALSE, algorithm = 3,resolution=resolution)
  return(Seu)
}



IGF131914 <- Norm(IGF131914)
IGF131914 <- RunTSNE(IGF131914,reduction = 'lsi')
DimPlot(RunTSNE(IGF131914,reduction='lsi',dims=2:50))
IGF131914 <- UMAP(IGF131914,dims=2:30,resolution = 0.2)
DimPlot(IGF131914,label=TRUE)+NoLegend()

IGF128020 <- Norm(IGF128020)
IGF128020 <- UMAP(IGF128020,dims=2:30)

IGF133741 <- Norm(IGF133741)
IGF133741 <- UMAP(IGF133741)
DimPlot(IGF133741,label=TRUE)+NoLegend()

IGF133740 <- Norm(IGF133740)
IGF133740 <- UMAP(IGF133740,resolution =0.7)
DimPlot(IGF133740,label=TRUE)+NoLegend()

geneAc <- GeneActivity(IGF133740)

markers0 <- FindMarkers(IGF133740,assay='activity',ident.1=0)
markers1 <- FindMarkers(IGF133740,assay='activity',ident.1=1)
markers2 <- FindMarkers(IGF133740,assay='activity',ident.1=2)
markers3 <- FindMarkers(IGF133740,assay='activity',ident.1=3)
markers4 <- FindMarkers(IGF133740,assay='activity',ident.1=4)

geneAc <- GeneActivity(IGF133741)
IGF133740[['activity']] <- CreateAssayObject(counts=geneAc)
IGF133741 <- NormalizeData(IGF133740,assay='activity')
FeaturePlot(object = IGF133740,features = c('ERG','CLDN5','VWF'),pt.size = 0.1)
FeaturePlot(object = IGF133740,features = c('AQP4','GFAP','SLC1A2','RFX4'),pt.size = 0.1)
FeaturePlot(object = IGF133740,features = c('NOTCH3','PDGFRB'),pt.size = 0.1)
FeaturePlot(object = IGF133740,features = c('MBP','ST18','SLC24A2','OLIG2','SOX10','OLIG1','MOBP'),pt.size = 0.1)
FeaturePlot(object = IGF133740,features = c('YTHDC1','TRA2B','RBFOX3','MAP2'),pt.size = 0.1)


#ASTROCYTES
FeaturePlot(object = IGF128020,features = c('SLC1A2','GFAP','SLC1A3','S100B','RFX4','ALDH1L1'),pt.size = 0.1)
#oligos
FeaturePlot(object = IGF128020,features = c('MBP','ST18','SLC24A2','OLIG2','SOX10','OLIG1','MOBP'),pt.size = 0.1)
#MICROGLIA
FeaturePlot(object = IGF128020,features = c('DOCK8','HS3ST4','IRF8'),pt.size = 0.1)
#NEURONS
FeaturePlot(object = IGF128020,features = c('YTHDC1','TRA2B'),pt.size = 0.1)
#ENDOTHELIAL
FeaturePlot(object = IGF128020,features = c('CLDN5'),pt.size = 0.1)
FeaturePlot(object = IGF128020,features = c('NOTCH3','RGS5','ABCC9'),pt.size = 0.1)
FeaturePlot(object = IGF128020,features = c('RBFOX3','MAP2'),pt.size = 0.1)

#CAPILLARY
FeaturePlot(object = IGF128020,features = c('FLT1','PECAM1','ERG','PODXL'),pt.size = 0.1)

markers0 <- FindMarkers(IGF128020,assay='activity',ident.1=0)
markers1 <- FindMarkers(IGF128020,assay='activity',ident.1=1)
markers2 <- FindMarkers(IGF128020,assay='activity',ident.1=2)
markers3 <- FindMarkers(IGF128020,assay='activity',ident.1=3)

geneAc2 <- GeneActivity(IGF131914)
IGF131914[['activity']] <- CreateAssayObject(counts=geneAc2)
IGF131914 <- NormalizeData(IGF131914,assay='activity')

markers131_0 <- FindMarkers(IGF131914,assay='activity',ident.1=0)
markers131_1 <- FindMarkers(IGF131914,assay='activity',ident.1=1)

IGF131914 <- RunTFIDF(IGF131914)
IGF131914 <- FindTopFeatures(IGF131914,min.cutoff = 'q0')
IGF131914 <- RunSVD(IGF131914)
DepthCor(IGF131914)
DimPlot(RunUMAP(IGF131914,reduction='lsi',dims=2:10,n.neighbors = 10))
DimPlot(IGF131914,label=TRUE)+NoLegend()
