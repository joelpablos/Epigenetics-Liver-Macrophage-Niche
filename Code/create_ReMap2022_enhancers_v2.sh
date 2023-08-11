#!/bin/sh

# Downsload mm39 UCSC genome assembly
#genomepy install mm39 -p UCSC --genomes_dir ../Data/Genomes/
#gzip ../Data/Genomes/mm39/mm39.fa

# Download annotation in HGNC gene symbols
wget ../Data/Genomes/mm39/https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_name_aliases&col=gd_pub_eg_id&col=gd_pubmed_ids&col=md_eg_id&col=md_refseq_id&col=md_ensembl_id&col=md_agr&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit

# Create directory to create the enhancers dataset for mm39
mkdir -p ../Data/REMAP2022/mm39

# Download all transcription factor ChIP-seq peaks for mm39 from the ReMap2022 databsase
wget -P ../Data/REMAP2022/mm39 https://remap.univ-amu.fr/storage/remap2022/mm39/MACS2/remap2022_all_macs2_mm39_v1_0.bed.gz
zcat ../Data/REMAP2022/mm39/remap2022_all_macs2_mm39_v1_0.bed.gz | head 

# Extract the summits of each peak
## If you want to remove the random regions and just keep the chromosomes add this command to the next line of code: grep -E '^chr[0-9XYM]+\b'
echo "Extracting the summits of each peak and extending them 25 bp up- and down-stream (and sorting them)..."
awk -v OFS='\t' '{$1;$2-=25;$3+=25; print}' <( zcat ../Data/REMAP2022/mm39/remap2022_all_macs2_mm39_v1_0.bed.gz | awk '{gsub(" +", "\t"); print}' | cut -f 1,7,8 ) > ../Data/REMAP2022/mm39/summit_extended.bed
sort -k1V -k2n -k3n ../Data/REMAP2022/mm39/summit_extended.bed > ../Data/REMAP2022/mm39/summit_extended_sort.bed

# Based on this file, generate a coverage bedGraph
## If you want to remove the random regions and just keep the chromosomes add this command to the next line of code: grep -E '^chr[0-9XYM]+\b'
echo "Generating coverage bedgraph based on bed file..."
bedtools genomecov -i ../Data/REMAP2022/mm39/summit_extended_sort.bed -g <(cat ../Data/Genomes/mm39/UCSC/mm39.fa.sizes | awk '{gsub(" +", "\t"); print}' | sort -k1V -k2n) -bg > ../Data/REMAP2022/mm39/enhancers.bedgraph
echo "Bedgraph generated succesfully. See here the head -n 50 ..."
head -n 50 ../Data/REMAP2022/mm39/enhancers.bedgraph

# Run peak calling using bdgpeakcall from MACS2 with l=50 and g=10
echo "Running MACS2 peak calling using bdgpeakcall with l=50 and g=10:"
## 1st: Run peak calling with cut-off c = 4
echo  "	- 1st: Running peak calling with cut-off c = 4 ..."
macs2 bdgpeakcall -i ../Data/REMAP2022/mm39/enhancers.bedgraph -l 50 -g 10 -c 4 -o ../Data/REMAP2022/mm39/enhancer_macs2_c4_peaks.bed
wc -l ../Data/REMAP2022/mm39/enhancer_macs2_c4_peaks.bed

## 2nd: Run peak calling again with cut-off c = 30
echo "	- 2nd: Running peak calling again with cut-off c = 30 ..."
macs2 bdgpeakcall -i ../Data/REMAP2022/mm39/enhancers.bedgraph -l 50 -g 10 -c 30 -o ../Data/REMAP2022/mm39/enhancer_macs2_c30_peaks.bed
wc -l ../Data/REMAP2022/mm39/enhancer_macs2_c30_peaks.bed

# Combine all peaks from c=30  with all peaks of c=4 that did not overlap with the peaks of c=30
bedtools subtract -a ../Data/REMAP2022/mm39/enhancer_macs2_c4_peaks.bed -b ../Data/REMAP2022/mm39/enhancer_macs2_c30_peaks.bed -A > ../Data/REMAP2022/mm39/c4_peaks_nonoverlapping.bed
cat ../Data/REMAP2022/mm39/enhancer_macs2_c30_peaks.bed ../Data/REMAP2022/mm39/c4_peaks_nonoverlapping.bed > ../Data/REMAP2022/mm39/peaks_macs2_overlap.bed
echo "All peaks from c=30 were combined with all peaks of c=4 that did not overlap with the peaks of c=30 succesfully. Length:"
wc -l ../Data/REMAP2022/mm39/peaks_macs2_overlap.bed

# Remove regions on chrM
grep -v "chrM" ../Data/REMAP2022/mm39/peaks_macs2_overlap.bed > ../Data/REMAP2022/mm39/peaks_macs2_overlap_filtered.bed
rm ../Data/REMAP2022/mm39/peaks_macs2_overlap.bed

# Extend summit of the peaks 100 bp up- and downstream to generate final enhancer collection of 1,268,775 putative enhancers of 200 bp
awk '{print $1"\t"$2+$10-100"\t"$2+$10+100}' ../Data/REMAP2022/mm39/peaks_macs2_overlap_filtered.bed | tail -n +2 | sort -k1V -k2n -k3n > ../Data/REMAP2022/mm39/ReMap2022_mm39_enhancers_v2.bed
echo "Enhancers regions of 200 bp long created!"
wc -l ../Data/REMAP2022/mm39/ReMap2022_mm39_enhancers_v2.bed
head ../Data/REMAP2022/mm39/ReMap2022_mm39_enhancers_v2.bed

