# ---------------------------------------------------------
# Integration pipeline for single-cell RNA-seq datasets
# Input: Preprocessed scRNA-seq Seurat objects (.rds)
# Output: One integrated Seurat object (.rds)
# ---------------------------------------------------------

rm(list=ls())
gc()
library('Seurat')
library('SeuratDisk')
library('dplyr')
library('patchwork')
library(ggplot2)
library('tidyverse')
library(SingleCellExperiment)

# --- Load datasets ---
# Read preprocessed datasets into a list object (E.list)
# Example:
# E1 <- readRDS("E1.rds")
# E2 <- readRDS("E2.rds")
# E.list <- list(E1, E2, ...)
# rm (E1, E2, ...)

# --- Integration setup ---
# Select integration features
features <- SelectIntegrationFeatures(object.list = E.list, nfeatures = 1e6)

# Load TE list (transposable element subfamilies)
TE_info <- read.table("mouse_te.txt", sep = "\t")
TE_info$V1 <- gsub("_", "-", TE_info$V1)
TE.list <- TE_info$V1


# Split anchor features into genes vs. TEs
anchor.features.TE <- intersect(features,TE.list) 
anchor.features <-SelectIntegrationFeatures(object.list = E.list,nfeatures=5000) 
anchor.features.gene <-setdiff(anchor.features,anchor.features.TE)[1:3000] # keep only the top 3000 genes


# --- Integration ---
E.anchors <- FindIntegrationAnchors(object.list = E.list,anchor.features = anchor.features)
rm(E.list); gc()
E.combine<-IntegrateData(anchorset = E.anchors)
rm(E.anchors);gc()

DefaultAssay(E.combine) <- "integrated"

# Scale and run PCA based on highly variable gene expression level
integraed_genes<-setdiff(row.names(E.combine),TE.list)
E.combine <- ScaleData(E.combine,vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(features =integraed_genes,assay = 'integrated') 

# --- Dimensional reduction, clustering and UMAP --- 
ElbowPlot(E.combine)
DimensionPicker<-function(seurat_obj){
  pct <- seurat_obj [["pca"]]@stdev / sum( seurat_obj [["pca"]]@stdev)*100 # Computes the proportion of variance explained by each principal component
  cumu <- cumsum(pct) #the cumulative sum of the percentages (how much of the total variance is explained as more principal components are included)
  co1 <- which(cumu > 90 & pct < 5)[1] #the indices of principal components where the cumulative variance explained is greater than 90% and the individual percentage of variance explained by the PC is less than 5%.
  co2 <- which(cumu > 95 & pct < 5)[1] #the indices of principal components where the cumulative variance explained is greater than 95% and the individual percentage of variance explained by the PC is less than 5%.
  co3 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 #Computes the difference in percentage variance explained between consecutive principal components & selects the largest drop
  output<-c(co1,co2,co3)
  return(output)
}
dim <- DimensionPicker(E.combine)[3]

E.combine<-FindNeighbors(E.combine,dims = 1:dim) %>% 
  FindClusters() %>% 
  RunUMAP(dims = 1:dim)


# --- Cell type annotation ---
DefaultAssay(E.combine)<-'RNA'
E.combine<-NormalizeData(E.combine) %>% ScaleData()
Idents(E.combine)<-'Time'

mks<-FindAllMarkers(E.combine,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
mks<-mks[mks$p_val_adj<0.05, ]
write.csv(mks,'gene_mks.csv')
# Assign cell type identity to each cluster in the metadata

# Visualization
DimPlot(E.combine, group.by = "celltype", label = TRUE, repel = TRUE) + 
  DimPlot(E.combine,group.by = 'Phase') + 
  DimPlot(E.combine, group.by = "Time")

# Save integrated object
output_path <- '.'
saveRDS(E.combine,file = file.path(output_path,"Integrated.rds"))

# --- Re-do dimensional reduction with TE expression only (within EN lineage cells) ---
E.enlineage <- subset(E.combine,
                      celltype %in% c("NE & RGC", "IPC", "EN"))

integraed_TEs<-setdiff(row.names(E.enlineage),TE.list)

E.enlineage <- ScaleData(E.enlineage) %>% 
  RunPCA(features =integraed_TEs,assay = 'integrated') 

ElbowPlot(E.enlineage)
dim <- DimensionPicker(E.enlineage)[3]

E.enlineage<-FindNeighbors(E.enlineage,dims = 1:dim) %>% FindClusters() %>% RunUMAP(dims = 1:dim)

DimPlot(E.enlineage, group.by = "celltype") + 
  DimPlot(E.enlineage,group.by = 'Phase') + 
  DimPlot(E.enlineage, group.by = "Time")

# Save integrated object
output_path <- '.'
saveRDS(E.combine,file = file.path(output_path,"Enlineage.rds"))


