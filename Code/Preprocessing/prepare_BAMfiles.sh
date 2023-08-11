#!/bin/sh

# Before performing ananse binding BAM files need to be sorted and indexed
# rm ../Data/ATAC/HSC/BAM/*tmp*
# rm ../Data/CNR/HSC/BAM/*tmp*

echo "Sorting and indexing BAM files..."

# cell_types=("KC" "HSC" "LSEC" "Hep" "Neutro")
cell_types=("Neutro")

for cell in "${cell_types[@]}"; do
	# Define the list of BAM files to be sorted and indexed
	bam_files=(
		"../Data/ATAC/${cell}/BAM_galaxy/${cell}_ATAC_REP1.bam"
		"../Data/ATAC/${cell}/BAM_galaxy/${cell}_ATAC_REP2.bam"
		"../Data/CNR/${cell}/BAM_galaxy/${cell}_H3K27Ac_REP1.bam"
		"../Data/CNR/${cell}/BAM_galaxy/${cell}_H3K27Ac_REP2.bam"
		)
	 bam_files=(
                "../Data/CNR/${cell}/BAM_galaxy/${cell}_H3K27Ac_REP1.bam"
                "../Data/CNR/${cell}/BAM_galaxy/${cell}_H3K27Ac_REP2.bam"
                )

	# Loop over each BAM file, sort it and index it
	for bam_file in "${bam_files[@]}"; do
		# Generate the corresponding sorted BAM file name
		echo "$(date) - ${bam_file} file sorting and indexing..."
	
		# Sort the BAM file
		bam_sorted="${bam_file%.*}_sort.bam"
		samtools sort "$bam_file" --threads 16 -o "$bam_sorted"
		echo "$(date) - ${bam_file} file sorted successfully!"

		# Index the sorted BAM file
		bam_indexed="${bam_sorted%.*}_index.bam"
		samtools index "$bam_sorted" -@ 16

		# Print a success message
		echo "$(date) - ${bam_indexed} file sorted and indexed successfully!"
	done
done
