# ---------------------------------------------------------
# Time-course clustering of TE subfamilies with TCseq
# Input: Seurat object (.rds) and TE subfamily list (mouse_te.txt)
# Output: Cluster trajectory plots (PDF)
# ---------------------------------------------------------

rm(list = ls())
gc()

library(TCseq)
library(dplyr)
library(Seurat)

# --- Load inputs ---
# Load TE list and seurat_obj
# seurat_obj <- readRDS("E.combine.rds")
# seurat_obj <- subset(seurat_obj, celltype=="NE & RGCs")
# TE_info <- read.table("mouse_te.txt", sep = "\t", header = FALSE)
# TE_info$V1 <- gsub("_", "-", TE_info$V1)  # standardize TE names to match matrix rows
# TE.list <- TE_info$V1

# --- Pseudobulk expression ---
DefaultAssay(seurat_obj)<-'RNA'
pb.expdata<-NormalizeData(seurat_obj) %>% 
  ScaleData(features = TE.list) %>% 
  AverageExpression(group.by = 'celltype',assays = 'RNA',features = TE.list) %>% 
  as.data.frame() 
colnames(pb.expdata)<-c('E10.5', 'E11.5','E12.5','E14.5')

# --- Time-course clustering ---
pbmatrix<-as.matrix(pb.expdata)

# Soft clustering with fuzzy cmeans
tca<-timeclust(pbmatrix,algo = 'cm',k=12,standardize = T) 

# Plot cluster trajectories
color_palette <- c("#F7FBFF","#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C")
plot<-timeclustplot(tca,
                    value = 'z-score',
                    membership.color = color_palette,
                    categories = "Timepoints", 
                    cols = 4)

# Save plots
plot_name<-'TCseq-cluster'
ggsave(paste0(plot_name,'.pdf'),plot,height = 210,width = 297,units = 'mm', dpi=96)

# --- Optional:Save objects ---
# saveRDS(tca, file = "TCseq_tca_timepoints.rds")

