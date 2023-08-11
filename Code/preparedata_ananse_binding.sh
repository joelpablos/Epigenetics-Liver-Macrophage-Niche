#!/bin/sh

# Before performing ananse binding BAM files need to be indexed
echo "Indexing BAM files..."
## ATAC-seq
samtools index -@ 8 ../Data/ATAC/KC/BAM/KC_ATAC_1.bam ../Data/ATAC/KC/BAM/KC_ATAC_1.bai
echo "	KC_ATAC_1.bam file indexed successfully "
samtools index -@ 8 ../Data/ATAC/KC/BAM/KC_ATAC_2.bam ../Data/ATAC/KC/BAM/KC_ATAC_2.bai
echo "	KC_ATAC_2.bam file indexed successfully "

## CUT&RUN
samtools index -@ 8 ../Data/CNR/KC/BAM/KC_H3K27Ac_1.bam ../Data/CNR/KC/BAM/KC_H3K27Ac_1.bai
echo "	KC_H3K27Ac_1.bam file indexed successfully "
samtools index -@ 8 ../Data/CNR/KC/BAM/KC_H3K27Ac_2.bam ../Data/CNR/KC/BAM/KC_H3K27Ac_2.bai
echo "	KC_H3K27Ac_2.bam file indexed successfully "
