# ---------------------------------------------------------
# hdWGCNA pipeline for integrated EN lineage cells
# Input: Integrated Seurat object (.rds)
# Output: Integrated Seurat object with hdWGCNA information(.rds)
# Reference: https://smorabit.github.io/hdWGCNA/index.html
# ---------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(UCell)

# --- Global settings ---
theme_set(theme_cowplot())   # Use cowplot theme
set.seed(12345)              # Reproducibility
allowWGCNAThreads(nThreads = 16) # Parallelization

# --- Input ---
input_path <- "."
seurat_obj <- readRDS(file.path(input_path, "Enlineage.rds"))
DefaultAssay(seurat_obj) <- "RNA"

# Restrict to TE features if TE.list is available
TE_info <- read.table("mouse_te.txt", sep = "\t")
TE_info$V1 <- gsub("_", "-", TE_info$V1)
TE.list <- TE_info$V1
seurat_obj <- seurat_obj[rownames(seurat_obj) %in% TE.list, ]

# --- Data preprocessing ---
seurat_obj<-NormalizeData(seurat_obj)%>%ScaleData()


# --- Setup for WGCNA ---
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  group.by='celltype',
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Co-expression" # the name of the hdWGCNA experiment
)

# --- Construct metacells ---
seurat_obj<- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c('celltype','Time'), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 100, # nearest-neighbors parameter
  max_shared = 20, # maximum number of shared cells between two metacells
  ident.group = 'celltype',# set the Idents of the metacell seurat object
)

seurat_obj<- NormalizeMetacells(seurat_obj)
gc()

# --- Expression setup ---
seurat_obj<-SetDatExpr(
  seurat_obj,
  group_name = c('NE & RGC', "IPC",'EN'),
  group.by='celltype',
  assay = 'RNA',
  slot = 'data'
)

# --- Soft-thresholding ---
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table,20)

# --- Co-expression network construction ---
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power= 8, 
  setDatExpr=FALSE,
  minModuleSize = 5,
  overwrite_tom = T,
  deepSplit = 4,
  tom_name = 'Celltype-coexpression' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='Celltype hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)

# --- Compute module eigengenes ---
seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj))
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Time",
  wgcna_name = "Co-expression" 
)

MEs <- GetMEs(seurat_obj, harmonized=FALSE)
hMEs<-GetMEs(seurat_obj)

# --- Compute connectivity and hub genes ---
seurat_obj <- ModuleConnectivity(
  seurat_obj
)
PlotKMEs(seurat_obj, ncol=5)

# get the module assignment table:
modules <- GetModules(seurat_obj)
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)

# --- Module scoring (with UCell) ---
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 'all',
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
wrap_plots(plot_list, ncol=3)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = T # depending on Seurat vs UCell for gene scoring
)
wrap_plots(plot_list, ncol=3)

# --- Correlation analysis ---
ModuleCorrelogram(seurat_obj,features = 'hMEs')


# Add hMEs to metadata for visualization
MEs <- GetMEs(seurat_obj, harmonized=T)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

DotPlot(seurat_obj, features=mods, group.by = 'celltype',dot.scale = 10)+
  coord_flip() +
  scale_color_gradientn(values = seq(0,1,0.2),colors = (c('#176BA0','black','#EE9A3A')))