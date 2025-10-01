# ---------------------------------------------------------
# Preprocessing pipeline for single-cell RNA-seq (Seurat)
# Input: scTE .h5ad file
# Output: Filtered, normalized, clustered Seurat object (.rds)
# ---------------------------------------------------------

rm(list=ls())
gc()
library('Seurat')
library('SeuratDisk')
library(gprofiler2)
library('dplyr')
library(ggplot2)
library(DropletUtils)
library("DoubletFinder")
library('tidyverse')
library(SingleCellExperiment)

# --- Load dataset ---
sample = "sample_name" # replace with the real sample file name
input_path <- '.' # replace with the real input directory
h5ad_file <- file.path(input_path,paste0(sample,"/",sample,".h5ad"))
# Convert the AnnData to Seurat data
Convert(h5ad_file, dest = "h5seurat")
seurat_obj<-LoadH5Seurat(file.path(input_path,paste0(sample,".h5seurat")))
seurat_obj<-CreateSeuratObject(seurat_obj@assays$RNA)

# Optional: if pre-saved RDS exists
# seurat_obj <- readRDS(file.path(sample, paste0(sample, ".rds")))


# --- QC Step 1: Remove empty droplets ---
rank<-barcodeRanks(seurat_obj@assays[["RNA"]]@counts)
uniq<-!duplicated(rank$rank)

plot(rank$rank[uniq], rank$total[uniq], log="xy",xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(rank)$inflection,col="darkgreen", lty=2)
abline(h=metadata(rank)$knee,col="dodgerblue", lty=2)
legend("topright", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=0.5)

cells<-rank@rownames[rank@listData[["total"]]>=17400] # replace the threshold with the knee value
seurat_obj<-subset(seurat_obj, cells = cells) 


# --- QC Step 2: Remove low-quality cells ---
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj, pattern="^mt-")
seurat_obj[["percent.rp"]]<-PercentageFeatureSet(seurat_obj, pattern="^Rp[sl]")
seurat_obj@meta.data$log10GenesPerUMI<-log10(seurat_obj$nFeature_RNA)/log10(seurat_obj$nCount_RNA)

VlnPlot(seurat_obj,features=c("nFeature_RNA","nCount_RNA",'log10GenesPerUMI','percent.mt','percent.rp'),ncol=5)

# Select thresholds based on violin plot visualization for each sample
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA < 8000 &
                       percent.mt < 3.5 &
                       percent.rp < 40 &
                       log10GenesPerUMI>0.8) 

# --- Cell cycle scoring ---
cc.genes.updated.2019 # Human cell cycle genes
# Map to mouse cell cycle genes
s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

# --- Normalization & PCA ---
seurat_obj<-NormalizeData(seurat_obj) %>% 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>% 
  SCTransform(variable.features.n = 60000,vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(features = rownames(seurat_obj))

# --- Dimensionality reduction & clustering ---
DimensionPicker<-function(seurat_obj){
  pct <- seurat_obj [["pca"]]@stdev / sum( seurat_obj [["pca"]]@stdev)*100 # Computes the proportion of variance explained by each principal component
  cumu <- cumsum(pct) #the cumulative sum of the percentages (how much of the total variance is explained as more principal components are included)
  co1 <- which(cumu > 90 & pct < 5)[1] #the indices of principal components where the cumulative variance explained is greater than 90% and the individual percentage of variance explained by the PC is less than 5%.
  co2 <- which(cumu > 95 & pct < 5)[1] #the indices of principal components where the cumulative variance explained is greater than 95% and the individual percentage of variance explained by the PC is less than 5%.
  co3 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 #Computes the difference in percentage variance explained between consecutive principal components & selects the largest drop
  output<-c(co1,co2,co3)
  return(output)
}

ElbowPlot(seurat_obj)
dim = DimensionPicker(seurat_obj)[3]
# Use the third output of DimensionPicker for FindNeighbors
seurat_obj<-FindNeighbors(seurat_obj,dims = 1:dim) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:dim)

# --- Doublet detection ---
DoubleRate=ncol(seurat_obj)*8*1e-6
annotations <- seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)          
nExp_poi <- round(DoubleRate*ncol(seurat_obj))  

seurat_obj<- doubletFinder_v3(seurat_obj, 
                              PCs = 1:30, pN = 0.25, pK = 0.01, 
                              nExp = nExp_poi, reuse.pANN = FALSE,sct = T)
seurat_obj<-subset(seurat_obj,DF.classifications_0.25_0.01_101!='Doublet') # adjust DF.classifications_xxx according to actual results

# --- Save results ---
# seurat_obj$Time <- "E10.5"  # optional metadata
output_path<-'.'
saveRDS(seurat_obj,file = file.path(output_path,paste0(sample, ".rds")))

