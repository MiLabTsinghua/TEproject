#!/bin/bash

# Enable immediate exit on error and define an error trap
set -e  # Stop script on error
trap 'echo "❌ Aborted due to an error in command: \"$BASH_COMMAND\"" && exit 1' ERR  # Print error message on failure

# Function to display a progress bar
progress_bar() {
    local duration=$1
    local step=0
    echo -n "["
    while [ $step -lt $duration ]; do
        echo -n "#"
        sleep 1
        ((step++))
    done
    echo "] Done!"
}

# Prompt user for required inputs
read -p "Enter project path (projPath): " projPath
read -p "Enter histName: " histName


echo "========== Step 1: Merge BAM Files =========="
# Check if input BAM files exist before merging
if [[ ! -f "${projPath}/alignment/bam/${histName}-rep1.sorted.bam" || ! -f "${projPath}/alignment/bam/${histName}-rep2.sorted.bam" ]]; then
    echo "❌ Error: One or both input BAM files are missing! Merging aborted."
    exit 1
fi

echo "✅ Input BAM files found. Proceeding with merging..."
samtools merge -@ 16 -o ${projPath}/alignment/bam/${histName}.bam \
    ${projPath}/alignment/bam/${histName}-rep1.sorted.bam \
    ${projPath}/alignment/bam/${histName}-rep2.sorted.bam 

# Index BAM file 
samtools index -@ 16 ${projPath}/alignment/bam/${histName}.bam
progress_bar 100  # Simulate progress

echo "========== Step 2: Generate Bigwig Files =========="
bamCoverage -b $projPath/alignment/bam/${histName}.bam \
-o $projPath/alignment/bigwig/${histName}.bw --normalizeUsing RPKM -p 8
progress_bar 100
read -p "Enter final name for bigwig file: " finalName
mv $projPath/alignment/bigwig/${histName}.bw $projPath/alignment/bigwig/$finalName.bw
echo "Bigwig file saved as $projPath/alignment/bigwig/$finalName.bw"

echo "========== 🎉 All steps completed successfully! 🎉 =========="