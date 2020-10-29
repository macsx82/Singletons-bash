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
    T)
        # stats concat step
        (echo "CHROM START END N_SITES TD";cat ${file_path}/${chr}.*.Tajima_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.TajimaD
        # we also need a BED formatted version for overlapping purposes (without header)
        (cat ${file_path}/${chr}.*.Tajima_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.TajimaD.bed
        
    ;;
    H)

        # concat Het stats
        (echo "CHROM START END TOT_HET N_SITES MEAN_HET MEAN_HET_RATIO";cat ${file_path}/${chr}.*.het_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.HET
        # we also need a BED formatted version for overlapping purposes (without header)
        (cat ${file_path}/${chr}.*.het_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.HET.bed

    ;;
    R)

        # STEP 4: ROH stats

    ;;
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
    I)
        # concat IHS stats: 
        (echo "CHROM START END N_SITES AVG_IHS AVG_SD_IHS";cat ${file_path}/${chr}.*.win_IHS_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.IHS
        # we also need a BED formatted version for overlapping purposes (without header)
        (cat ${file_path}/${chr}.*.win_IHS_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.IHS.bed

    ;;
    C)
        # concat CADD
        (echo "CHROM START END BIN_0_10 BIN_10_20 BIN_20 N_CADD_ANN N_TOT_WIN";cat ${file_path}/${chr}.*.CADD_win_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.CADD
        # we also need a BED formatted version for overlapping purposes (without header)
        (cat ${file_path}/${chr}.*.CADD_win_stats | sort -g -k3,3 -k2,2)| tr " " "\t" > ${out_path}/${chr}.ALL.CADD.bed

    ;;
esac
done
