#!/bin/bash
# ---------------------------------------------------------
# CUT&Tag Processing Pipeline
# Input: Raw FASTQ files (paired-end)
# Output: BAM, BED, BedGraph, BigWig (RPKM normalized)
# Dependencies: Trim Galore, Bowtie2, Samtools, Bedtools, deeptools
# ---------------------------------------------------------

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

ref="$HOME/ref/mm10_ref/mm10"
spikeInRef="$HOME/ref/Ecoli/ecoli"
chromSize="$HOME/ref/mm10_ref/mm10.chrom.sizes.txt"
cores=16

echo "========== Step 1: Trimming Reads =========="
# Trim reads
while true; do
    read -p "Run Trimgalore? [y/n/other] " trim
    if [ "$trim" == "y" ]; then

        echo "Running Trim Galore with adapter type: $adapter ..."
        trim_galore --phred33 --fastqc --stringency 10 --gzip --length 50 --max_n 10 -q 30 \
        --nextera -j 8 --paired \
        ${projPath}/Rawdata/${histName}/${histName}_R1.fq.gz \
        ${projPath}/Rawdata/${histName}/${histName}_R2.fq.gz \
        --output_dir ${projPath}/Rawdata/${histName}/

        progress_bar 100  # Simulate progress
        fq_file1=${projPath}/Rawdata/${histName}/${histName}_R1_val_1.fq.gz
        fq_file2=${projPath}/Rawdata/${histName}/${histName}_R2_val_2.fq.gz
        break
    elif [ "$trim" == "n" ]; then
        echo "Skipping Trim Galore..."
        fq_file1=${projPath}/Rawdata/${histName}/${histName}_R1.fq.gz
        fq_file2=${projPath}/Rawdata/${histName}/${histName}_R2.fq.gz
        break
    else
        echo "Invalid input. Please enter 'y' or 'n'."
    fi
done


echo "========== Step 2: Alignment to Reference Genome =========="
mkdir -p ${projPath}/alignment/sam/bowtie2_summary ${projPath}/alignment/bam ${projPath}/alignment/bed ${projPath}/alignment/bedgraph
echo "Aligning reads to reference genome..."
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} \
-x ${ref} -1 ${fq_file1} -2 ${fq_file2} | samtools view -@ ${cores} -bS - > ${projPath}/alignment/sam/${histName}_bowtie2.bam \
2>> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
progress_bar 100


echo "========== Step 3: Alignment to Spike-in Genome =========="
echo "Aligning reads to spike-in genome (E. coli)..."
bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant \
--phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ${fq_file1} -2 ${fq_file2} | \
samtools view -@ ${cores} -bS - > $projPath/alignment/sam/${histName}_bowtie2_spikeIn.bam \
2>> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
progress_bar 100


echo "========== Step 4: Calculating Sequencing Depth =========="
seqDepthDouble=$(samtools view -@ 16 -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.bam | wc -l)
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth
progress_bar 100


echo "========== Step 5: Fragment Size Distribution =========="
mkdir -p $projPath/alignment/sam/fragmentLen
samtools view -@ 16 -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.bam | \
awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' \
>$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt
progress_bar 100


echo "========== Step 6: Filtering and Bed Conversion =========="
samtools view -@ 16 -b -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.bam \
>$projPath/alignment/bam/${histName}_bowtie2.mapped.bam
bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.mapped.bam -bedpe \
>$projPath/alignment/bed/${histName}_bowtie2.bed
progress_bar 100


echo "========== Step 7: Filtering by Fragment Length and Chromosome =========="
awk '$1==$4 && $6-$2 < 1000 {print $0}' \
$projPath/alignment/bed/${histName}_bowtie2.bed \
>$projPath/alignment/bed/${histName}_bowtie2.clean.bed
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.fragments.bed
progress_bar 100

echo "========== Step 8: Spike-in Calibration =========="
if [[ "$seqDepth" -gt "1" ]]; then
    mkdir -p $projPath/alignment/bedgraph
    scale_factor=$(echo "10000 / $seqDepth" | bc -l)
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor \
    -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed \
    -g $chromSize \
    >$projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
    progress_bar 100
fi


echo "========== Step 9: Generating BigWig Files =========="
mkdir -p $projPath/alignment/bigwig
samtools sort -@ 16 -o $projPath/alignment/bam/${histName}.sorted.bam \
$projPath/alignment/bam/${histName}_bowtie2.mapped.bam
samtools index -@ 16 $projPath/alignment/bam/${histName}.sorted.bam
bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam \
-o $projPath/alignment/bigwig/${histName}_RPKM.bw --normalizeUsing RPKM -p ${cores}
progress_bar 100

echo "========== Step 10: File Cleaning =========="
read -p "Do you want to remove intermediate files? [y/n]: " clean
if [[ "clean" == "y" ]]; then
    rm $projPath/alignment/sam/${histName}_bowtie2.bam
    rm $projPath/alignment/sam/${histName}_bowtie2_spikeIn.bam
    if [ "$trim" == "y" ]; then
        rm fq_file1 fq_file2 
    fi
    echo "Intermediate files removed successfully!"
else
    echo "Intermediate files kept."
fi
progress_bar 100


echo "========== 🎉 All steps completed successfully! 🎉 =========="