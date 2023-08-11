#!/bin/sh

echo "$(date) - Performing multi-pairwise bidirectional ANANSE influence over all cells ..."

database="" # write here the name of the database used e.g. ("gimme39")
db_prefix="" # db_prefix="/${database}_"

cell_types1=("KC_KG" "Mono_KG")
cell_types2=("HSC" "LSEC" "Hep" "KC" "Mono_KG" "KC_KG")

# Loop over ANANSE pipeline for the choosen cells
for source in "${cell_types1[@]}"; do
	echo "$(date) - ANANSE influence: ${source} to all cells: "
	mkdir ../Results/ANANSE_influence${database}/
	mkdir ../Results/ANANSE_influence${database}/${source}/
	mkdir ../Results/ANANSE_plot${database}/
	mkdir ../Results/ANANSE_plot${database}/${source}/

	for target in "${cell_types2[@]}"; do
		source_cell=$(echo "${source}" | sed 's/_.*//g')
		target_cell=$(echo "${target}" | sed 's/_.*//g')

		if [ "$source_cell" != "$target_cell" ]; then			
			if true;then 
			echo "	$(date) - performing ANANSE influence from ${source} to ${target} cell ..."
			ananse influence  -s ../Results/ANANSE_network/${source}/${db_prefix}ANANSE_network.tsv \
					-t ../Results/ANANSE_network/${target}/${db_prefix}ANANSE_network.tsv \
					-d ../Data/ATLAS/MouseStSt/specificScore/PRESTO/VS_DEGs/PRESTO_DEGs_${source_cell}vs${target_cell}.tsv \
					-o ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence.tsv \
					-f \
					 --select-after-join \
					-i 1_000_000 \
					-n 16
			echo "  $(date) - ANANSE influence: from ${source} to ${target} cell performed successfully!!"
			
			fi
			echo "  $(date) - plotting ANANSE diff.network from ${source} to ${target} cell ..."
			ananse plot ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence.tsv \
				--diff-network ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence_diffnetwork.tsv \
				--n-tfs 30 \
				#--node-placement twopi \
				-t png \
				-o ../Results/ANANSE_plot${database}/${source}/${source}vs${target}
			echo "  $(date) - ANANSE plot from ${source} to ${target} cell generated successfully!!"
			
		fi
	done
done


genes=("Spic" "Rbpj" "Nr1h3" "Rxra" "Ets1" "Tcf7l2" "Irf7")
# cell_types=("KC" "LSEC")
echo -e "\n$(date) - ANANASE influence statistics for pairwise [${cell_types[@]}]: "

for source in "${cell_types[@]}"; do
	echo -e "---------------------------------------------------------------------------------------------------------------------------\n ${source} pairwise combinations:"
	for target in "${cell_types[@]}"; do
                if [ "$source" != "$target" ]; then
			echo -e "\n - Number of edges AND influential TFs: ${source}vs${target} AND ${target}vs${source}"
			wc -l ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence*
			wc -l ../Results/ANANSE_influence${database}/${target}/${target}vs${source}_influence*
			
			echo -e "\t - Top 10 most influential TFs"
			echo -e "\t\t - ${source}vs${target}:"
			cat ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence.tsv | cut -f 1,2 | sort -k 2n | tail -n 10 | tac
			echo -e "\t\t - ${target}vs${source}:"
			cat ../Results/ANANSE_influence${database}/${target}/${target}vs${source}_influence.tsv | cut -f 1,2 | sort -k 2n | tail -n 10 | tac
			# echo -e "\t\t - Common:"
			# comm <(cat ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence.tsv | cut -f 1,2 | sort -k 2n | tail -n 10 | cut -f 1 | sort) \
			#      <(cat ../Results/ANANSE_influence${database}/${target}/${target}vs${source}_influence.tsv | cut -f 1,2 | sort -k 2n | tail -n 10 | cut -f 1 | sort) \
			     
			if false; then
			echo -e "\n\t - Looking for specific genes [${genes[@]}] in ${source}vs${target} AND ${target}vs${source}"
			for gene in "${genes[@]}"; do
				echo -e "\t----------------------------------------------------------------------\n\t\t - Results for ${gene}:"
				if false; then
					echo -e "\n\t\t\t - ${source}vs${target} diff.network:"
					grep -i ${gene} ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence_diffnetwork.tsv | cut -f 1,2,3 | sort -k 3n | tail -n 10 | tac
					echo -e "\t\t\t - ${target}vs${source} diff.network:"
					grep -i ${gene} ../Results/ANANSE_influence${database}/${target}/${target}vs${source}_influence_diffnetwork.tsv | cut -f 1,2,3 | sort -k 3n | tail -n 10 | tac
					echo -e "\t\t\t - ${source}vs${target} AND ${target}vs${source} diff.network (comm):"
					comm <(grep -i ${gene} ../Results/ANANSE_influence${database}/${source}/${source}vs${target}_influence_diffnetwork.tsv | cut -f 1,2,3 | sort -k 3n | tail -n 10 | cut -f 1,2 | sort) <(grep -i ${gene} ../Results/ANANSE_influence/${database}/${target}/${target}vs${source}_influence_diffnetwork.tsv | cut -f 1,2,3 | sort -k 3n | tail -n 10 | cut -f 1,2 | sort)
				fi
				if false; then
					echo -e "\n\t\t\t - ${source}vs${target} influence TFs:"
                                	echo -e "\t\t\t\t$(grep -i ${gene} ../Results/ANANSE_influence/${database}/${source}/${source}vs${target}_influence.tsv | cut -f 1,2)"
                                	echo -e "\t\t\t - ${target}vs${source} influence TFs:"
                                	echo -e "\t\t\t\t$(grep -i ${gene} ../Results/ANANSE_influence/${database}/${target}/${target}vs${source}_influence.tsv | cut -f 1,2)"
				fi
			done
			fi
		fi
	done
done

