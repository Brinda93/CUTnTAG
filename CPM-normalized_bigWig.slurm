#!/bin/bash -l
#SBATCH --job-name=Sort_CPM_normalize
#SBATCH --output=logs/Sort_CPM_normalize.%A_%a.out
#SBATCH --error=logs/Sort_CPM_normalize.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --account=lab_dovat
#SBATCH --array=1-20

set -euo pipefail

# -------------------------
# Directories
# -------------------------
BASE_OUT="old_script1"
MAPPED_HG="${BASE_OUT}/Mapped_Hg"
SORTED_DIR="${BASE_OUT}/Sorted_Bam"
NORM_DIR="${BASE_OUT}/Normalized_CPM"
LOGS_DIR="${BASE_OUT}/logs"

mkdir -p "${SORTED_DIR}" "${NORM_DIR}" "${LOGS_DIR}"

# -------------------------
# Threads & genome
# -------------------------
THREADS=8
HG38_CHROM_SIZES="/gpfs/Labs/Dovat/bxp5423/genome/Genome/Hg38/hg38.chrom.sizes"

# Load required modules
module load samtools
module load deeptools/3.5.1

# -------------------------
# Detect all BAM files
# -------------------------
BAMS=($(ls ${MAPPED_HG}/*_trimmed.bam | xargs -n1 basename | sed 's/_trimmed.bam//' | sort))
NUM_SAMPLES=${#BAMS[@]}

echo "Detected $NUM_SAMPLES BAM files"

if [[ $SLURM_ARRAY_TASK_ID -gt $NUM_SAMPLES ]]; then
    echo "Array index $SLURM_ARRAY_TASK_ID exceeds number of samples ($NUM_SAMPLES)"
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID-1))
SAMPLE=${BAMS[$INDEX]}

INPUT_BAM="${MAPPED_HG}/${SAMPLE}_trimmed.bam"
SORTED_BAM="${SORTED_DIR}/${SAMPLE}_sorted.bam"
OUTPUT_BIGWIG="${NORM_DIR}/${SAMPLE}_CPM.bw"

echo "Processing sample: $SAMPLE"
echo "Input BAM: $INPUT_BAM"

# -------------------------
# Step 1: Sort BAM file by coordinate
# -------------------------
echo "Sorting BAM file..."
samtools sort \
    -@ ${THREADS} \
    -m 3G \
    -o "${SORTED_BAM}" \
    "${INPUT_BAM}" \
    2> "${LOGS_DIR}/${SAMPLE}_sort.log"

echo "Sorted BAM created: ${SORTED_BAM}"

# -------------------------
# Step 2: Index sorted BAM file
# -------------------------
echo "Indexing sorted BAM file..."
samtools index -@ ${THREADS} "${SORTED_BAM}"

# -------------------------
# Step 3: Count total reads
# -------------------------
echo "Counting total reads..."
TOTAL_READS=$(samtools view -@ ${THREADS} -c -F 260 "${SORTED_BAM}")
echo "Total reads: ${TOTAL_READS}"

# Calculate scaling factor (CPM = reads per million)
SCALE_FACTOR=$(echo "scale=10; 1000000 / ${TOTAL_READS}" | bc)
echo "CPM scaling factor: ${SCALE_FACTOR}"

# Save statistics
echo -e "Sample\tTotal_Reads\tScaling_Factor" > "${LOGS_DIR}/${SAMPLE}_CPM_stats.txt"
echo -e "${SAMPLE}\t${TOTAL_READS}\t${SCALE_FACTOR}" >> "${LOGS_DIR}/${SAMPLE}_CPM_stats.txt"

# -------------------------
# Step 4: Generate CPM-normalized bigWig using deepTools
# -------------------------
echo "Generating CPM-normalized bigWig..."
bamCoverage \
    --bam "${SORTED_BAM}" \
    --outFileName "${OUTPUT_BIGWIG}" \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2913022398 \
    --numberOfProcessors ${THREADS} \
    --extendReads \
    2> "${LOGS_DIR}/${SAMPLE}_bamCoverage.log"

echo "✅ Sorting and CPM normalization completed for ${SAMPLE}"
echo "Sorted BAM: ${SORTED_BAM}"
echo "CPM bigWig: ${OUTPUT_BIGWIG}"
