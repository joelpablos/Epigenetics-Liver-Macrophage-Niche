#!/bin/sh

ananse network  ../Results/ANANSE_binding/binding.h5 \
                 -e ../Data/ATLAS/MouseStSt/KC_pseudobulk.tsv \
                 -c tpm
                 -g ../Data/Genomes/mm39 \
                 -o ../Results/ANANSE_network/ANANSE_network.tsv

echo "ANANSE network has been created successfully!!"
