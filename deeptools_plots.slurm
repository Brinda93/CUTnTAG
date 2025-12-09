module load deeptools/3.5.1
module load bedtools/2.31.0

 

cat CT4_CKDL220002265-1a_HK55MDSX3_L2_peaks.narrowPeak CT5_CKDL220002266-1a_HK55MDSX3_L2_peaks.narrowPeak CT6_CKDL220002267-1a_HK55MDSX3_L2_peaks.narrowPeak | sort -k1,1 -k2,2n  | bedtools merge -i - > WT_E2_C4_merged_peaks.bed


computeMatrix scale-regions     -R WT_K8_K13_merged_peaks.bed    -S ../Normalized_CPM/CT4_CKDL220002265-1a_HK55MDSX3_L2_CPM.bw ../Normalized_CPM/CT5_CKDL220002266-1a_HK55MDSX3_L2_CPM.bw ../Normalized_CPM/CT6_CKDL220002267-1a_HK55MDSX3_L2_CPM.bw -b 3000 -a 3000 --regionBodyLength 5000     --skipZeros     --numberOfProcessors 8     -o matrix_genes_TSS_TES_motl4.gz

plotHeatmap -m matrix_genes_TSS_TES_motl4.gz    -out heatmap_TSS_to_TES_Molt4.png     --colorMap RdYlBu_r     --samplesLabel "WT" "C2" "E4"     --regionsLabel "Genes"     --plotTitle "Molt4"     --startLabel "TSS"     --endLabel "TES" 


plotProfile -m matrix_genes_TSS_TES_motl4.gz  -out profile_WT_C2_E4.png     --plotTitle "Molt4"     --samplesLabel "WT" "C2" "E4"     --colors blue green red     --plotType se     --perGroup 

or


#!/bin/bash
#SBATCH --job-name=deeptools_peaks
#SBATCH --output=deeptools_peaks_%j.out
#SBATCH --error=deeptools_peaks_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=compute

module load deeptools/3.5.1
module load bedtools/2.31.0

cd /gpfs/Labs/Dovat/Rloop_Jun2022/Rloop_Peaks_Feb2022/CUT_Tag_Jan_2022_CS/raw_files/result_original_script_for_cross_check/old_script1

# Step 1: Convert peaks to genes with peaks
# Find genes that have peaks in their promoter regions
bedtools closest -a Peaks/CT1_CKDL220002262-1a_HHM3TDSX3_L3_peaks.narrowPeak \
    -b /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed -d | awk '$NF < 3000' | cut -f1-3 | \
    bedtools intersect -a /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed -b - -u > genes_with_WT_peaks.bed

bedtools closest -a Peaks/CT2_CKDL220002263-1a_HHM3TDSX3_L3_peaks.narrowPeak \
    -b /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed -d | awk '$NF < 3000' | cut -f1-3 | \
    bedtools intersect -a /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed -b - -u > genes_with_K8_peaks.bed

bedtools closest -a Peaks/CT3_CKDL220002264-1a_HK55MDSX3_L2_peaks.narrowPeak \
    -b /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed -d | awk '$NF < 3000' | cut -f1-3 | \
    bedtools intersect -a /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed -b - -u > genes_with_K13_peaks.bed

echo "Gene filtering completed!"

# Step 2: Compute matrix for genes WITH peaks
computeMatrix scale-regions \
    -S Normalized_CPM/CT1_CKDL220002262-1a_HHM3TDSX3_L3_CPM.bw \
       Normalized_CPM/CT2_CKDL220002263-1a_HHM3TDSX3_L3_CPM.bw \
       Normalized_CPM/CT3_CKDL220002264-1a_HK55MDSX3_L2_CPM.bw \
    -R genes_with_WT_peaks.bed \
       genes_with_K8_peaks.bed \
       genes_with_K13_peaks.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --skipZeros \
    -p 8 \
    -o matrix_genes_with_peaks_nalm6.gz \
    --samplesLabel WT K8 K13

echo "computeMatrix completed successfully!"

# Step 3: Plot heatmap
plotHeatmap -m matrix_genes_with_peaks_nalm6.gz \
    -out heatmap_genes_with_peaks_nalm6.png \
    --samplesLabel WT K8 K13 \
    --regionsLabel "Genes_with_WT_peaks" "Genes_with_K8_peaks" "Genes_with_K13_peaks" \
    --plotTitle "Nalm6 - Genes with Peaks" \
    --colorMap RdYlBu_r \
    --startLabel "TSS" \
    --endLabel "TES" \
    --xAxisLabel "gene distance (bp)" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --zMin 0

echo "plotHeatmap completed successfully!"

# Print summary statistics
echo "=== Peak-to-Gene Summary ==="
echo "Genes with WT peaks: $(wc -l < genes_with_WT_peaks.bed)"
echo "Genes with K8 peaks: $(wc -l < genes_with_K8_peaks.bed)"
echo "Genes with K13 peaks: $(wc -l < genes_with_K13_peaks.bed)"
