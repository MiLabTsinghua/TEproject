#!/bin/bash
# ---------------------------------------------------------
# Mapping Pipeline for Single-Cell Transposable Element Expression Data
# Input: Raw FASTQ files (paired-end) from single-cell RNA-seq
# Output: AnnData object with TE and gene expression matrices
# Dependencies: Cell Ranger, scTE
# ---------------------------------------------------------

# Enable immediate exit on error and define an error trap
set -e  # Stop script on error
trap 'echo "❌ Aborted due to an error in command: \"$BASH_COMMAND\"" && exit 1' ERR  # Print error message on failure
# Prompt user for required inputs
read -p "Enter project path: " projPath
read -p "Enter sample name: " sampleName

# 1. Cell Ranger count
# Note: Update the transcriptome and fastqs paths according to your setup
echo "Running Cell Ranger count for sample: $sampleName ..."
cellranger count --id=${sampleName} \
 --transcriptome=$HOME/ref/mm10_ref/cellranger_mm10 \
 --fastqs=${projPath}/${sampleName}/ \
 --sample=${sampleName} \
 --create-bam=true \
 --include-introns=true

# 2. Run scTE to generate TE and gene expression matrices
# Note: Update the genome and TE annotation paths according to your setup
echo "Running scTE for sample: $sampleName ..."
scTE -i ${sampleName}/outs/possorted_genome_bam.bam \
 -o ${projPath}/${sampleName}_scTE \
 -x $HOME/scTE/mm10.exclusive.idx  \
 -CB CB -UMI UB \
 --hdf5 True -p 8

echo "Mapping and TE quantification completed for sample: $sampleName"
echo "Output files are located in: ${projPath}/${sampleName}_scTE"



