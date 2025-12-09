#!/bin/bash -l
#SBATCH --job-name=Rloop_ALL
#SBATCH --output=old_script/logs/Rloop_ALL.%A_%a.out
#SBATCH --error=old_script/logs/Rloop_ALL.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=compute
#SBATCH --account=lab_dovat
#SBATCH --array=1-20   # <-- update after checking number of samples

set -euo pipefail

# -------------------------
# Directories
# -------------------------
INPUT_DIR="/gpfs/Labs/Dovat/Rloop_Jun2022/Rloop_Peaks_Feb2022/CUT_Tag_Jan_2022_CS/raw_files"
BASE_OUT="old_script1"
MAPPED_HG="${BASE_OUT}/Mapped_Hg"
TAG_DIR="${BASE_OUT}/TagDir"
PEAKS_DIR="${BASE_OUT}/Peaks"
LOGS_DIR="${BASE_OUT}/logs"
mkdir -p "${MAPPED_HG}" "${TAG_DIR}" "${PEAKS_DIR}" "${LOGS_DIR}"

# -------------------------
# Threads & genome
# -------------------------
THREADS=8
HG38_INDEX="/gpfs/Labs/Dovat/bxp5423/genome/Genome/Hg38/hg38.bowtie_index"
HG38_CHROM_SIZES="/gpfs/Labs/Dovat/bxp5423/genome/Genome/Hg38/hg38.chrom.sizes"

# Load required modules
module load bowtie2/2.4.2
module load samtools
module load homer/4.10
export PATH=/ri/shared/modules7/UCSCTools/406/bin/:$PATH
module load macs2/2.2.9.1

# -------------------------
# Detect all samples
# -------------------------
SAMPLES=($(ls ${INPUT_DIR}/*_1.fq.gz | sed 's/_1.fq.gz//' | xargs -n1 basename | sort))
NUM_SAMPLES=${#SAMPLES[@]}
echo "Detected $NUM_SAMPLES samples"

if [[ $SLURM_ARRAY_TASK_ID -gt $NUM_SAMPLES ]]; then
    echo "Array index $SLURM_ARRAY_TASK_ID exceeds number of samples ($NUM_SAMPLES)"
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID-1))
SAMPLE=${SAMPLES[$INDEX]}
READS_1="${INPUT_DIR}/${SAMPLE}_1.fq.gz"
READS_2="${INPUT_DIR}/${SAMPLE}_2.fq.gz"

echo "Processing sample: $SAMPLE"
echo "Reads: $READS_1 and $READS_2"

# -------------------------
# Verify tools
# -------------------------
for tool in bowtie2 makeTagDirectory makeUCSCfile macs2 samtools; do
    command -v $tool >/dev/null 2>&1 || { echo "Error: $tool not found in PATH"; exit 1; }
done

# -------------------------
# Bowtie2 alignment
# -------------------------
ALN_SAM="${MAPPED_HG}/${SAMPLE}_trimmed.sam"
bowtie2 -x "${HG38_INDEX}" \
    --trim3 100 \
    -p "${THREADS}" \
    -1 "${READS_1}" \
    -2 "${READS_2}" \
    -S "$ALN_SAM" \
    2> "${LOGS_DIR}/${SAMPLE}_hg38_trimmed_alignment.log"

ln -sf "$ALN_SAM" "${MAPPED_HG}/${SAMPLE}.sam"

# -------------------------
# Convert SAM -> BAM
# -------------------------
ALN_BAM="${MAPPED_HG}/${SAMPLE}_trimmed.bam"
samtools view -@ "${THREADS}" -bS "$ALN_SAM" > "$ALN_BAM"

# -------------------------
# HOMER tag directory & bigWig
# -------------------------
makeTagDirectory "${TAG_DIR}/${SAMPLE}_hg38/" "$ALN_BAM"
makeUCSCfile "${TAG_DIR}/${SAMPLE}_hg38/" \
    -bigWig "$HG38_CHROM_SIZES" \
    -o "${TAG_DIR}/${SAMPLE}_hg38/${SAMPLE}_hg38"

# -------------------------
# MACS2 peak calling
# -------------------------
macs2 callpeak \
    -t "$ALN_BAM" \
    -f BAMPE \
    -q 0.05 \
    -n "$SAMPLE" \
    --outdir "${PEAKS_DIR}/" \
    --nomodel \
    2> "${LOGS_DIR}/${SAMPLE}_macs2.log"

# -------------------------
# Generate summary
# -------------------------
SUMMARY_FILE="${LOGS_DIR}/${SAMPLE}_summary.txt"
{
echo "=== Alignment Summary ==="
echo "Sample: $SAMPLE"
echo "Date: $(date)"
echo "Bowtie2 version: $(bowtie2 --version | head -n1)"
echo ""
echo "Total alignments in BAM:"
samtools view -c "$ALN_BAM"
echo ""
echo "Bowtie2 log:"
cat "${LOGS_DIR}/${SAMPLE}_hg38_trimmed_alignment.log"
} > "$SUMMARY_FILE"

echo "✅ Completed $SAMPLE"

