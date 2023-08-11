#!/bin/sh

echo "$(date) - Performing gimme scan over ReMap2022 mm39 enhancer landscape..."
# gimme scan -Tz --gc -g -N 20 ../Data/Genomes/mm39 ../Data/REMAP2022/mm39/ReMap2022_mm39_enhancers_v2.bed > ../Data/REMAP2022/mm39/SCAN.tsv
echo "$(date) - gimme scan over enhancers performed successfully!!"

# cell_types=("KC" "HSC" "LSEC" "Hep")
# cell_types=("KC_KG" "Mono_KG")
CNRorChIPseq=("CNR")

# Loop over ANANSE pipeline for the choosen cells
for cell in "${cell_types[@]}"; do

	# ananse binding with gimme motif2factor mm39 database adapted from gimme.vertebrates.v5
	ananse binding -A ../Data/nf-core/atacseq/${cell}/bwa/merged_library/${cell}_REP1.mLb.clN.sorted.bam \
			  ../Data/nf-core/atacseq/${cell}/bwa/merged_library/${cell}_REP2.mLb.clN.sorted.bam \
	               -H ../Data/${CNRorChIPseq}/${cell}/BAM_galaxy/${cell}_H3K27Ac_REP1.bam \
			  ../Data/${CNRorChIPseq}/${cell}/BAM_galaxy/${cell}_H3K27Ac_REP2.bam \
	               -g ../Data/Genomes/mm39 \
	               -r ../Data/REMAP2022/mm39/ReMap2022_mm39_enhancers_v2.bed.gz \
	               -p ../Data/GimmeMotifs/gimmevertebratev5_mm39/mm39.gimme.vertebrate.v5.0.pfm \
	               --pfmscorefile ../Data/REMAP2022/mm39/SCAN.tsv \
	               -n 20 \
		       -o ../Results/ANANSE_binding/${cell}/

	echo "$(date) - ANANSE binding.h5 has been created successfully!!"

	ananse network  ../Results/ANANSE_binding/${cell}/binding.h5 \
			-e ../Data/ATLAS/MouseStSt/${cell}_pseudobulk.tsv \
			-c cpm \
			-g ../Data/Genomes/mm39 \
			--full-output \
			-n 20 \
			-o ../Results/ANANSE_network/${cell}/ANANSE_network.tsv
	
	echo "$(date) - ANANSE network has been created successfully!!"
	
	# We make the network compatible with Cytoscape visualization
	cat ../Results/ANANSE_network/${cell}/ANANSE_network.tsv | awk -F'\t' 'NR==1{split($1,a,"_"); print a[1]"\t"a[2]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6} NR>1{split($1,a,"—"); print a[1]"\t"a[2]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ../Results/ANANSE_network/${cell}/ANANSE_network_CYTOSCAPE.tsv

	echo "$(date) - ANANSE network for CYTOSCAPE generated successfully!! ../Results/ANANSE_network/${cell}/ANANSE_network_CYTOSCAPE.tsv"
	
	echo "$(date) - ANANSE view: here you have the binding probabilities of the first 10 regions and motifs"
	ananse view ../Results/ANANSE_binding/${cell}/binding.h5 -n 10

	echo "$(date) - ANANSE view: ${cell} regions, TFs and binding probabilities will be extracted and saved in ../Results/ANANSE_view/${cell}/"
	echo "... (loading) ..."
	ananse view ../Results/ANANSE_binding/${cell}/binding.h5 -lr -o ../Results/ANANSE_view/${cell}/regions_binding.txt 
	echo "$(date) - ANANSE view: List of regions in the ../Results/ANANSE_binding/${cell}/binding.h5 file saved: ../Results/ANANSE_view/${cell}/regions_binding.txt"

	ananse view ../Results/ANANSE_binding/${cell}/binding.h5 -lt -o ../Results/ANANSE_view/${cell}/TFs_binding.txt
        echo "$(date) - ANANSE view: List of transcription factors in the ../Results/ANANSE_binding/${cell}/binding.h5 file saved: ../Results/ANANSE_view/${cell}/TFs_binding.txt"

	ananse view ../Results/ANANSE_binding/${cell}/binding.h5 -a -o ../Results/ANANSE_view/${cell}/TFactivity_binding.txt
        echo "$(date) - ANANSE view: List of transcription factors activity scores in the ../Results/ANANSE_binding/${cell}/binding.h5 file saved: ../Results/ANANSE_view/${cell}/TFactivity_binding.txt"

	ananse view ../Results/ANANSE_binding/${cell}/binding.h5 -o ../Results/ANANSE_view/${cell}/regionsVsTFs_binding.tsv
	echo "$(date) - ANANSE view: binding probabilities for all factors in wide format (one column per TF) have been extracted: ../Results/ANANSE_view/${cell}/regionsVsTFs_binding.tsv"

done

exit

source="KC"
target="Neutro"

echo "$(date) - ANANSE influence: performing ananse influence between source (${source}) and target (${target}) network"
ananse influence  -s ../Results/ANANSE_network/${source}/ANANSE_network.tsv \
                    -t ../Results/ANANSE_network/${target}/ANANSE_network.tsv \
                    -d ../Data/ATLAS/MouseStSt/DESeq2_${source}vs${target}.tsv \
		    -n 12 \
                    -o ../Results/ANANSE_influence/${source}vs${target}influence.tsv

echo "$(date) - ANANSE influence: ../Results/ANANSE_influence/${source}vsHSC_influence.tsv and ../Results/ANANSE_influence/KCvsHSC_influence_diffnetwork.tsv generated successfully!!"

echo "$(date) - ANANSE plot: plotting top differential TFs between ${source} vs ${target}"
ananse plot ../Results/ANANSE_influence/${source}vs${target}_influence.tsv \
                --diff-network ../Results/ANANSE_influence/${source}vs${target}_influence_diffnetwork.tsv \
                --n-tfs 30 \
                -t jpeg \
                -f \
                -o ../Results/ANANSE_plot/
echo "$(date) - ANANSE plot: plot of top 30 differential TFS between ${source} and ${target} generated succesfully!!"

echo "$(date) - ANANSE workflow finished! Thanks for using me! (By Joel Pablos Martin - IRC-VIB Gent, Belgium)"

