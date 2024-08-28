library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(readr)
library(ggplot2)
library(dplyr)
library(harmony)

set.seed(23)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
end(annotations) <- (start(annotations) + 500)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"


peaks125165 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF125165/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks128020 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF128020/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks131914 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF131914/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks133738 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133737/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks133740 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133740/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks133741 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133741/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks133742 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133742/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))

peaks133743 <- makeGRangesFromDataFrame(read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133743/outs/peaks.bed",
  col.names = c("chr", "start", "end")
))


combinedPeaks <- reduce(x = c(peaks125165,peaks128020,peaks131914,peaks133738,peaks133740,peaks133741,peaks133742,peaks133743))
peakwidths <- width(combinedPeaks)
combinedPeaks <- combinedPeaks[peakwidths  < 10000 & peakwidths > 20]

md125165 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF125165/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md128020 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF128020/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md131914 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF131914/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md133738 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133738/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md133740 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133740/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md133741 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133741/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md133742 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133742/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md133743 <- read.table(
  file = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133743/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md125165 <- md125165[md125165$is__cell_barcode == 1, ]
md128020 <- md128020[md128020$is__cell_barcode == 1, ]
md131914 <- md131914[md131914$is__cell_barcode == 1, ]
md133738 <- md133738[md133738$is__cell_barcode == 1, ]
md133740 <- md133740[md133740$is__cell_barcode == 1, ]
md133741 <- md133741[md133741$is__cell_barcode == 1, ]
md133742 <- md133742[md133742$is__cell_barcode == 1, ]
md133743 <- md133743[md133743$is__cell_barcode == 1, ]

frags125165 <-  CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF125165/outs/fragments.tsv.gz",
  cells = rownames(md125165)
)

frags128020 <-  CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF128020/outs/fragments.tsv.gz",
  cells = rownames(md128020)
)

frags131914 <-  CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF131914/outs/fragments.tsv.gz",
  cells = rownames(md131914)
)

frags133738 <- CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133738/outs/fragments.tsv.gz",
  cells = rownames(md133738)
)

frags133740 <- CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133740/outs/fragments.tsv.gz",
  cells = rownames(md133740)
)

frags133741 <- CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133741/outs/fragments.tsv.gz",
  cells = rownames(md133741)
)

frags133742 <- CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133742/outs/fragments.tsv.gz",
  cells = rownames(md133742)
)

frags133743 <- CreateFragmentObject(
  path = "/rds/general/user/cg2723/home/SC/submitfiles/IGF133743/outs/fragments.tsv.gz",
  cells = rownames(md133743)
)


counts125165 <- FeatureMatrix(
  fragments = frags125165,
  features = combinedPeaks,
  cells = rownames(md125165)
)

counts128020 <- FeatureMatrix(
  fragments = frags128020,
  features = combinedPeaks,
  cells = rownames(md128020)
)

counts131914 <- FeatureMatrix(
  fragments = frags131914,
  features = combinedPeaks,
  cells = rownames(md131914)
)


counts133738 <- FeatureMatrix(
  fragments = frags133738,
  features = combinedPeaks,
  cells = rownames(md133738)
)

counts133740 <- FeatureMatrix(
  fragments = frags133740,
  features = combinedPeaks,
  cells = rownames(md133740)
)

counts133741 <- FeatureMatrix(
  fragments = frags133741,
  features = combinedPeaks,
  cells = rownames(md133741)
)

counts133742 <- FeatureMatrix(
  fragments = frags133742,
  features = combinedPeaks,
  cells = rownames(md133742)
)

counts133743 <- FeatureMatrix(
  fragments = frags133743,
  features = combinedPeaks,
  cells = rownames(md133743)
)

assay125165 <- CreateChromatinAssay(counts125165,fragments = frags125165)
IGF125165 <- CreateSeuratObject(assay125165, assay = "CnT", meta.data = md125165)
IGF125165$sample <- 'IGF125165'
IGF125165$Sample <- '1'
IGF125165$patient <- 'P70'
IGF125165 <- RenameCells(IGF125165,add.cell.id = 'IGF125165')

assay128020 <- CreateChromatinAssay(counts128020,fragments = frags128020)
IGF128020 <- CreateSeuratObject(assay128020, assay = "CnT", meta.data = md128020)
IGF128020$sample <- 'IGF128020'
IGF128020$Sample <- '2'
IGF128020$patient <- 'P70'
IGF128020 <- RenameCells(IGF128020,add.cell.id = 'IGF128020')

assay131914 <- CreateChromatinAssay(counts131914,fragments = frags131914)
IGF131914 <- CreateSeuratObject(assay131914, assay = "CnT", meta.data = md131914)
IGF131914$sample <- 'IGF131914'
IGF131914$Sample <- '3'
IGF131914$patient <- 'P70'
IGF131914 <- RenameCells(IGF131914,add.cell.id = 'IGF131914')

assay133738 <- CreateChromatinAssay(counts133738,fragments = frags133738)
IGF133738 <- CreateSeuratObject(assay133738, assay = "CnT", meta.data = md133738)
IGF133738$sample <- 'IGF133738'
IGF133738$Sample <- '4'
IGF133738$patient <- 'P83'
IGF133738 <- RenameCells(IGF133738,add.cell.id = 'IGF133738')

assay133740 <- CreateChromatinAssay(counts133740,fragments = frags133740)
IGF133740 <- CreateSeuratObject(assay133740, assay = "CnT", meta.data = md133740)
IGF133740$sample <- 'IGF133740'
IGF133740$Sample <- '5'
IGF133740$patient <- 'P83'
IGF133740 <-RenameCells(IGF133740,add.cell.id = 'IGF133740')

assay133741 <- CreateChromatinAssay(counts133741,fragments = frags133741)
IGF133741 <- CreateSeuratObject(assay133741, assay = "CnT", meta.data = md133741)
IGF133741$sample <- 'IGF133741'
IGF133741$Sample <- '6'
IGF133741$patient <- 'P83'
IGF133741 <- RenameCells(IGF133741,add.cell.id = 'IGF133741')

assay133742 <- CreateChromatinAssay(counts133742,fragments = frags133742)
IGF133742 <- CreateSeuratObject(assay133742, assay = "CnT", meta.data = md133742)
IGF133742$sample <- 'IGF133742'
IGF133742$Sample <- '7'
IGF133742$patient <- 'P94'
IGF133742 <-RenameCells(IGF133742,add.cell.id = 'IGF133742')

assay133743 <- CreateChromatinAssay(counts133743,fragments = frags133743)
IGF133743 <- CreateSeuratObject(assay133743, assay = "CnT", meta.data = md133743)
IGF133743$sample <- 'IGF133743'
IGF133743$Sample <- '8'
IGF133743$patient <- 'P94'
IGF133743 <-RenameCells(IGF133743,add.cell.id = 'IGF133743')


combined <- merge(
  x = IGF125165,
  y = c(IGF128020,IGF131914,IGF133738,IGF133740,IGF133741,IGF133742,IGF133743)
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined)+NoLegend()
ggsave('/rds/general/user/cg2723/home/SC/plotsnew/UMAP_nointegration.png', dpi = 1000, width = 8, height = 5)
DimPlot(combined, group.by = "Sample")
ggsave('/rds/general/user/cg2723/home/SC/plotsnew/UMAP_nointegration_bysample.png', dpi = 1000, width = 8, height = 5)

DimPlot(combined, group.by = "patient")
ggsave('/rds/general/user/cg2723/home/SC/plotsnew/UMAP_nointegration_bypatient.png', dpi = 1000, width = 8, height = 5)

Norm <- function(Seu){
  Seu <- RunTFIDF(Seu)
  Seu <- FindTopFeatures(Seu, min.cutoff = 'q0')
  Seu <- RunSVD(Seu)
  return(Seu)
}
IGF125165 <- Norm(IGF125165)
IGF128020 <- Norm(IGF128020)
IGF131914 <- Norm(IGF131914)
IGF133738 <- Norm(IGF133738)
IGF133740 <- Norm(IGF133740)
IGF133741 <- Norm(IGF133741)
IGF133742 <- Norm(IGF133742)

integration.anchors <- FindIntegrationAnchors(
  object.list = list(IGF125165,IGF128020,IGF131914,IGF133738,IGF133740,IGF133741,IGF133742),
  reduction = "rlsi",
  dims = 2:30,
  anchor.features = rownames(combined),
  k.anchor=10
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)
DimPlot(integrated)

hm.integrated <- RunHarmony(object = combined, group.by.vars = 'patient', reduction = 'lsi', assay.use = 'CnT', theta = 4, lambda = 2, project.dim = FALSE)
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated, group.by = 'patient', pt.size = 0.1)
DimPlot(hm.integrated, group.by = 'sample', pt.size = 0.1)
hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'harmony', dims = 2:30)
hm.integrated <- FindClusters(object = hm.integrated, verbose = FALSE, algorithm = 3)
markers <- FindMarkers(hm.integrated, assay = 'activity', ident.1 = )
DimPlot(hm.integrated)

PatientsByCluster <- data.frame(cluster = hm.integrated@active.ident,patient =hm.integrated@meta.data$patient)

PatientsByCluster %>% group_by(cluster,patient) %>% summarise(n=n())

ggplot(PatientsByCluster, aes(cluster,fill = patient ))+geom_bar(position='fill')

Annotation(IGF133741) <- annotations
geneAc <- GeneActivity(IGF133741)

Annotation(combined) <- annotations
geneAc <- GeneActivity(combined)

Annotation(hm.integrated) <- annotations
geneAc <- GeneActivity(hm.integrated)
hm.integrated[['activity']] <- CreateAssayObject(counts=geneAc)
DimPlot(combined, group.by = 'patient', pt.size = 0.1)
combined[['activity']] <- CreateAssayObject(counts=geneAc)
combined <- NormalizeData(combined,assay='activity')
FeaturePlot(object = combined,features = c('ERG','CLDN5','VWF','OCLN','ESAM'),pt.size = 0.1)
FeaturePlot(object = combined,features = c('AQP4','GFAP','SLC1A2','RFX4'),pt.size = 0.1)
FeaturePlot(object = combined,features = c('MBP','ST18','SLC24A2','OLIG2','SOX10','OLIG1','MOBP'),pt.size = 0.1)
FeaturePlot(object = combined,features = c('YTHDC1','TRA2B','RBFOX3','MAP2'),pt.size = 0.1)
FeaturePlot(object = combined,features = c('NOTCH3','PDGFRB','RGS5','CSPG4','TAGLN'),pt.size = 0.1)


FeaturePlot(object = hm.integrated,features = c('ERG','CLDN5','VWF','OCLN','ESAM'),pt.size = 0.1)
FeaturePlot(object = hm.integrated,features = c('ERG'),pt.size = 0.1, split.by = 'patient')
FeaturePlot(object = hm.integrated,features = c('AQP4','GFAP','SLC1A2','RFX4'),pt.size = 0.1)
FeaturePlot(object = hm.integrated,features = c('NOTCH3','PDGFRB','RGS5','CSPG4','TAGLN'),pt.size = 0.1)
FeaturePlot(object = hm.integrated,features = c('MBP','ST18','SLC24A2','OLIG2','SOX10','OLIG1','MOBP'),pt.size = 0.1)
FeaturePlot(object = hm.integrated,features = c('SPI1','AIF1','TMEM119','IRF8'),pt.size = 0.1)
FeaturePlot(object = hm.integrated,features = c('CDC45'),pt.size = 0.1)


DimPlot(object = hm.integrated,pt.size = 0.1, label = T)
markers11 <- FindMarkers(hm.integrated, assay = 'activity', ident.1 = 11)
