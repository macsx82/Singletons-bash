#!/usr/bin/env bash

###############################################
#23/04/2018
#
#Component for Pipeline script for the variants prioritization project
# general script used to concat back chunked vcf files or other files
#v 0.1
###############################################

#ARGS
file_path=$2
chr=$3
out_path=$4

#for each option we need to concat back the stats files generated

while getopts ":THRSIC" opt; do
    case ${opt} in
    S)
        # concat Singletons output
        (echo -e "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV";cat ${file_path}/${chr}.*-*.singletons | sort -g -k2,2| fgrep -v "ALLELE")| tr " " "\t" > ${out_path}/${chr}.ALL.samples.singletons
        # concat Singletons stats
        (echo "CHROM START END N_SING N_SITES SING_RATE W_size SING_RATE_by_W_size";cat ${file_path}/${chr}.*.singletons_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.singletons
        # we also need a BED formatted version for overlapping purposes (without header)
        (cat ${file_path}/${chr}.*.singletons_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.singletons.bed

        #add a cleaning and compressing step to reduce footprint of the generated data
        tar -czf ${file_path}/${chr}.singletons_stats.tar.gz ${file_path}/${chr}.*.singletons_stats
        tar -czf ${file_path}/${chr}.singletons.tar.gz ${file_path}/${chr}.*-*.singletons

        #clean single region data
        rm ${file_path}/${chr}.*.singletons_stats
        rm ${file_path}/${chr}.*-*.singletons
    ;;
esac
done
