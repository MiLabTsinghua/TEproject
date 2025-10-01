#!/bin/bash
# ---------------------------------------------------------
# Automated MACS2 peak calling and MSPC analysis for CUT&Tag
# Input: BAM files (sorted or name-sorted) from CUT&Tag mapping of one timepoint
# Output: MACS2 peak files and MSPC combined peaks
# Dependencies: MACS2, MSPC, Samtools
# ---------------------------------------------------------

# Enable immediate exit on error and define an error trap
set -e  # Stop script on error
trap 'echo "❌ Aborted due to an error in command: \"$BASH_COMMAND\"" && exit 1' ERR  # Print error message on failure

# Get timepoint for peaks: 
read -p "Enter timepoint: " timepoint

# Define histone modification types and patterns
declare -A MODS
MODS=( ["H3K4me3"]="H3K4me3" ["H3K4me1"]="H3K4me1" ["H3K27ac"]="H3K27ac" ["H3K9me3"]="H3K9me3" ["H3K27me3"]="H3K27me3" ["IgG"]="IgG")

# Initialize MOD_BAMS, only process original .sorted.bam and -rep*.bam files, and only add MOD if not empty
declare -A MOD_BAMS
for BAM in *.sorted.bam *.bam; do
    # 只处理支持的命名格式，且不是mapped.bam
    [[ "$BAM" == *.mapped.bam ]] && continue
    MOD=""
    if [[ "$BAM" =~ ^.*H3[^_]*_[0-9]+\.sorted\.bam$ ]]; then
        MOD=$(echo "$BAM" | awk -F'_' '{print $(NF-2)}')
    elif [[ "$BAM" =~ ^.*-H3[^-]*-rep[0-9]+\.bam$ ]]; then
        MOD=$(echo "$BAM" | awk -F'-' '{for(i=1;i<=NF;i++) if($i ~ /^H3/) print $i}')
    fi
    [[ -z "$MOD" ]] && continue
    MOD_BAMS[$MOD]+="$BAM "
done

# 0. Perform name sort for reps without .mapped.bam (skip if already exists)
for BAM in *.sorted.bam *.bam; do
    [[ "$BAM" == *.mapped.bam ]] && continue
    MOD=""
    NAME_SORTED_BAM=""
    if [[ "$BAM" =~ ^.*H3[^_]*_[0-9]+\.sorted\.bam$ ]]; then
        MOD=$(echo "$BAM" | awk -F'_' '{print $(NF-2)}')
        NAME_SORTED_BAM="${BAM/.sorted.bam/.mapped.bam}"
    elif [[ "$BAM" =~ ^.*-H3[^-]*-rep[0-9]+\.bam$ ]]; then
        MOD=$(echo "$BAM" | awk -F'-' '{for(i=1;i<=NF;i++) if($i ~ /^H3/) print $i}')
        NAME_SORTED_BAM="${BAM/.bam/.mapped.bam}"
    else
        echo "[Info] $BAM does not match supported naming pattern, skip."
        continue
    fi
    # Skip if .mapped.bam already exists
    if [ -f "$NAME_SORTED_BAM" ]; then
        echo "$NAME_SORTED_BAM already exists, skip this rep."
        continue
    fi
    echo "Sorting $BAM by read name..."
    samtools sort -n -o "$NAME_SORTED_BAM" "$BAM" -@ 16
done

# 1. MACS2 peak calling (ask each time whether to use IgG as control)
for BAM in *.mapped.bam; do
    # Skip IgG samples
    if [[ "$BAM" == *IgG* ]]; then
        continue
    fi
    # Process only supported naming patterns (.mapped.bam)
    MOD=""
    SAMPLE=""
    OUTDIR=""
    IGG_BAM=""
    IGG_NAME_SORTED_BAM=""
    if [[ "$BAM" =~ ^.*H3[^_]*_[0-9]+\.mapped\.bam$ ]]; then
        MOD=$(echo "$BAM" | awk -F'_' '{print $(NF-2)}')
        SAMPLE=$(basename "$BAM" .mapped.bam)
        OUTDIR="macs2_${SAMPLE}"
        IGG_BAM=$(echo "$BAM" | sed "s/${MOD}/IgG/")
        IGG_NAME_SORTED_BAM="${IGG_BAM/.mapped.bam/.mapped.bam}"
    elif [[ "$BAM" =~ ^.*-H3[^-]*-rep[0-9]+\.mapped\.bam$ ]]; then
        MOD=$(echo "$BAM" | awk -F'-' '{for(i=1;i<=NF;i++) if($i ~ /^H3/) print $i}')
        SAMPLE=$(basename "$BAM" .mapped.bam)
        OUTDIR="macs2_${SAMPLE}"
        IGG_BAM=$(echo "$BAM" | sed "s/${MOD}/IgG/")
        IGG_NAME_SORTED_BAM="${IGG_BAM/.mapped.bam/.mapped.bam}"
    else
        echo "[Info] $BAM does not match supported naming pattern, skip."
        continue
    fi
    mkdir -p $OUTDIR
    USE_CTRL="n"
    if [ -f "$IGG_BAM" ]; then
        read -p "Use $IGG_BAM as control for $BAM? [y/n]: " USE_CTRL
    fi

    # Skip IgG itself
    if [[ "$BAM" == *IgG* ]]; then
        continue
    fi
    SAMPLE=$(basename "$BAM" .sorted.bam)

    MOD=""
    if [[ "$BAM" =~ ^.*H3[^_]*_[0-9]+\.sorted\.bam$ ]]; then
        # 下划线分隔，修饰名为倒数第二段
        MOD=$(echo "$BAM" | awk -F'_' '{print $(NF-2)}')
    elif [[ "$BAM" =~ ^.*-H3[^-]*-rep[0-9]+\.sorted\.bam$ ]]; then
        # 中划线分隔，修饰名为带H3的段
        MOD=$(echo "$BAM" | awk -F'-' '{for(i=1;i<=NF;i++) if($i ~ /^H3/) print $i}')
    else
        echo "[Info] $BAM does not match supported naming pattern, skip."
        continue
    fi
    OUTDIR="macs2_${SAMPLE}"
    mkdir -p $OUTDIR
    # Auto-match IgG control (same batch same rep)
    IGG_BAM=$(echo "$BAM" | sed "s/${MOD}/IgG/")
    NAME_SORTED_BAM="${BAM/.sorted.bam/.name_sorted.bam}"
    IGG_NAME_SORTED_BAM="${IGG_BAM/.sorted.bam/.name_sorted.bam}"
    USE_CTRL="n"
    if [ -f "$IGG_BAM" ]; then
        read -p "Use $IGG_BAM as control for $BAM? [y/n]: " USE_CTRL
    fi
    # Set MACS2 parameters based on modification type
    MACS2_ARGS=""
    case "$MOD" in
        H3K4me3)
            MACS2_ARGS="--nomodel --extsize 147 -q 0.05 "
            ;;
        H3K4me1)
            MACS2_ARGS="--broad --broad-cutoff 0.1 "
            ;;
        H3K27ac)
            MACS2_ARGS="--nomodel --extsize 147 -q 0.05"
            ;;
        H3K9me3)
            MACS2_ARGS="--broad --broad-cutoff 0.1"
            ;;
        H3K27me3)
            MACS2_ARGS="--broad --broad-cutoff 0.1  "
            ;;
        *)
            MACS2_ARGS=""
            ;;
    esac
    # 组装命令（始终用-BAMPE）
    if [ "$USE_CTRL" == "y" ] && [ -f "$IGG_NAME_SORTED_BAM" ]; then
        macs2 callpeak -t $NAME_SORTED_BAM -c $IGG_NAME_SORTED_BAM -n $SAMPLE -f BAMPE --nolambda -g mm $MACS2_ARGS --outdir $OUTDIR
    else
        macs2 callpeak -t $NAME_SORTED_BAM -n $SAMPLE -f BAMPE --nolambda -g mm $MACS2_ARGS --outdir $OUTDIR
    fi
done

# 2. mspc
# For each modification, automatically match two replicates
for MOD in "${!MODS[@]}"; do
    # 获取该修饰的所有样本名
    SAMPLES=()
    for BAM in ${MOD_BAMS[$MOD]}; do
        SAMPLES+=( $(basename "$BAM" .sorted.bam) )
    done
    if [ ${#SAMPLES[@]} -lt 2 ]; then
        echo "[Warning] $MOD less than 2 replicates, skip IDR."
        continue
    fi
    # 取前两个rep
    REP1=${SAMPLES[0]}
    REP2=${SAMPLES[1]}
    if [ "${MODS[$MOD]}" == "narrow" ]; then
        TYPE="narrowPeak"
        PEAK1="macs2_${REP1}/${REP1}_peaks.narrowPeak"
        PEAK2="macs2_${REP2}/${REP2}_peaks.narrowPeak"
    else
        TYPE="broadPeak"
        PEAK1="macs2_${REP1}/${REP1}_peaks.broadPeak"
        PEAK2="macs2_${REP2}/${REP2}_peaks.broadPeak"
    fi
    OUT_IDR="${timepoint}_${MOD}_mspc"
    mspc --i $PEAK1 $PEAK2 -r tec -s 1e-4 -w 0.05 -o $OUT_IDR
done

