#!/bin/sh

# Generate the motif database for mm39 taking as reference the  non-redundant, 
# clustered database of known vertebrate motifs: gimme.vertebrate.v5.0
gimme motif2factors --new-reference mm39 --outdir ../Data/GimmeMotifs/mm39 --tmpdir ../tmp/ --genomes_dir ../Data/Genomes/
# or HOCOMOCOv11_MOUSE 
gimme motif2factors --new-reference mm39 --database HOCOMOCOv11_MOUSE --outdir ../Data/GimmeMotifs/HOMOCOCO_mm39 --tmpdir ../tmp/ --genomes_dir ../Data/Genomes/
