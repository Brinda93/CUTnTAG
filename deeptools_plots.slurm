module load deeptools/3.5.1
module load bedtools/2.31.0

computeMatrix reference-point -S CT1_CKDL220002262-1a_HHM3TDSX3_L3_CPM.bw CT2_CKDL220002263-1a_HHM3TDSX3_L3_CPM.bw CT3_CKDL220002264-1a_HK55MDSX3_L2_CPM.bw 
-R /gpfs/Labs/Dovat/BLISS_Oct2025_Ranjit/bliss102025/output/breaks/hotspot_analysis/genes.bed --beforeRegionStartLength 3000     --regionBodyLength 5000     --afterRegionStartLength 3000    
--missingDataAsZero     --samplesLabel "Nalm6_WT" "Nalm6-Ik-KO-K8" "Nalm6-IK-KO-K13"     -o /gpfs/Labs/Dovat/Rloop_Jun2022/Rloop_Peaks_Feb2022/Plots/matrix_Nalm6_genebody.gz     -p 8


plotHeatmap     -m matrix_Molt4_genedody.gz -out heatmap_Molt4_g
enebodies.png     --colorMap Reds     --sortRegions descend


plotProfile -m matrix_Nalm6_genebody.gz -out profile_Nalm6_genebody.png --perGroup --samplesLabel "Nalm6 WT" "Nalm6-Ik-KO-K8" "Nalm6-IK-KO-K13" --regionsLabel "Genes" 
--plotTitle "Nalm6 R-loop signal across gene bodies" --colors red blue green --averageType mean --plotType se --yAxisLabel "R-loop signal"
