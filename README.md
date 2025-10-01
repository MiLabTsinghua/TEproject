# TEproject

This repository contains scripts and documentation for analyzing transposable element (TE) regulation and expression using CUT&Tag and single-cell RNA-seq data, corresponding to the paper in preparation: [Title and link to be added].

## Project Overview

This project provides a complete analytical workflow for investigating the epigenetic regulation and expression patterns of TEs during neocortical development. The pipeline consists of the following major components:

- **CUT&Tag Data Processing**: From raw FASTQ files to chromatin accessibility analysis
- **Peak Calling**: Histone modification peak detection using MACS2
- **Single-cell RNA-seq Analysis**: Preprocessing, integration, and TE expression analysis
- **Co-expression Network Analysis**: TE co-expression module identification using hdWGCNA

## File Structure

```
TEproject/
├── README.md
├── bash/
|   ├── scTEMapping.sh        # Mapping of single-cell TE and gene expression data
│   ├── CUTTagMapping.sh      # CUT&Tag data alignment and processing
│   ├── MergeRep.sh           # Biological replicate merging
│   └── PeakCalling.sh        # Peak calling and IDR analysis
└── R/
|   ├── 00_Preprocessing.R    # Single-cell data preprocessing
|   ├── 01_Integration.R      # Multi-sample integration analysis
|   ├── 02_hdWGCNA.R         # Co-expression network analysis with hdWGCNA
|   └── 03_TCseq.R          # TE expression dynamics analysis with TCseq (to be added)
└── data/
    └── mouse_te.txt          # Mouse transposable element annotations
```

## Dependencies

### Bash Script Dependencies

- **Single-cell Processing**: Cell Ranger, scTE
- **Sequence Processing**: Trim Galore, FastQC
- **Alignment Tools**: Bowtie2, Samtools
- **Genomic Analysis**: Bedtools, deepTools
- **Peak Calling**: MACS2, MSPC

### R Script Dependencies

```r
# Core analysis packages
Seurat, SeuratDisk, dplyr, ggplot2, tidyverse

# Single-cell specific analysis
DropletUtils, DoubletFinder, SingleCellExperiment

# Functional annotation
gprofiler2

# Network analysis
WGCNA, hdWGCNA, UCell

# Visualization
cowplot, patchwork
```

## Citations

Please cite our paper once available.

## License

This project follows the MIT License.

## Contact

For questions or suggestions, please contact through GitHub Issues.