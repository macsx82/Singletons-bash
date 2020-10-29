#!/usr/bin/env bash


###############################################
#23/04/2018
#
#Pipeline script for the variants prioritization project
#v 0.1
###############################################

# Pipeline workflow:

# TajimaD stats
# Het stats
# ROH stats
# Singleton density
# IHS
# CADD score distribution   

# Each by chromosome
# Each can run in parallel

#define a function to generate bed files with region boundaries
#need to find the updated coordinates for each chr, since the absolute dist
base_bash_scripts="/home/cocca/scripts/bash_scripts"
######################
# ARGS:
# $1 is reserved for pipeline options
#chr number (TODO: parse chr name in vcf to check if the chr name is just the number or chrN)
chr=$2
#size in bp
win_size=$3 
out_folder=$4
pop_vcf=$5
m=$6
q=$7
#we need to account for a sample list,and for a bed file containing genes 
#The list will be an INCLUSION list, with a single column of samples's ids to include in the analyses
sample_list=$8
#this bed file is useless, most of the times, we need to check if it's still usefull to have this option
genes_bed=$9
# pipe_step=$7
# STEP 1: create a bed file for each chromosome to work with regions in parallel
# a) retrieve updated chromosome coordinates
chrom=(`mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B --skip-column-names -e "select chrom, 0, size as coords  from hg19.chromInfo where chrom = 'chr${chr}';"`)
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B --skip-column-names -e "select chrom, 0, size as coords  from hg19.chromInfo where chrom NOT LIKE 'chr___%' and chrom NOT LIKE 'chrUn_%' order by chrom;"

start=${chrom[1]}
end=${chrom[2]}

#chunk_outfile=${outfolder}/${pop}_${chr}_10.chunks
chunk_mode="range"

chunk_outfile="${out_folder}/${chr}_chunk_file.txt"
gene_stats_folder="${out_folder}/genes"

mkdir -p ${out_folder}
mkdir -p ${gene_stats_folder}


# echo "${@}"
while getopts ":THRSIC" opt; do
    case ${opt} in
    T)
        # echo "/home/cocca/scripts/bash_scripts/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}"
        ${base_bash_scripts}/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}

        #Now we'll use the chunk file to submit job arrays to calculate all the needed stats
        # STEP 2: calculate tajimaD stats
        tajD_out="${out_folder}/TajimaD"
        tajD_logs="${out_folder}/LOGS/TajimaD"

        mkdir -p ${tajD_out}
        mkdir -p ${tajD_logs}


        a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_TD.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${tajD_out}"|qsub -t 1-${a_size} -o ${tajD_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${tajD_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N tajD_chr${chr} -l h_vmem=${m} -q ${q}
   
        # stats concat step
        echo "${base_bash_scripts}/VPP_GR2018_pipeline_CONCAT.sh -T ${tajD_out} ${chr} ${tajD_out}"|qsub -o ${tajD_logs}/chr${chr}_\$JOB_ID_concat.log -e ${tajD_logs}/chr${chr}_\$JOB_ID_concat.e -V -N tajD_concat_chr${chr} -hold_jid tajD_chr${chr} -l h_vmem=2G -q ${q}

        if [[ ${genes_bed} != "" ]]; then
        #intersect the provided genes list with our stats
        echo "intersect the provided genes list with our stats..."
        echo "bedtools intersect -a ${genes_bed} -b ${tajD_out}/${chr}.ALL.TajimaD.bed -wo > ${gene_stats_folder}/${chr}.genes_TajimaD.bed" | qsub -o ${tajD_logs}/chr${chr}_\$JOB_ID_bed_intersect.log -e ${tajD_logs}/chr${chr}_\$JOB_ID_bed_intersect.e -V -N tajD_bed_intersect_chr${chr} -hold_jid tajD_concat_chr${chr} -l h_vmem=2G -q ${q}
        fi


    ;;
    H)
        # echo "/home/cocca/scripts/bash_scripts/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}"
        ${base_bash_scripts}/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}

        # STEP 3: Het stats
        het_out="${out_folder}/HET"
        het_logs="${out_folder}/LOGS/HET"

        mkdir -p ${het_out}
        mkdir -p ${het_logs}

        a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_HET.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${het_out}"|qsub -t 1-${a_size} -o ${het_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${het_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N het_chr${chr} -l h_vmem=${m} -q ${q}
        
        # stats concat step
        echo "${base_bash_scripts}/VPP_GR2018_pipeline_CONCAT.sh -H ${het_out} ${chr} ${het_out}"|qsub -o ${het_logs}/chr${chr}_\$JOB_ID_concat.log -e ${het_logs}/chr${chr}_\$JOB_ID_concat.e -V -N het_concat_chr${chr} -hold_jid het_chr${chr} -l h_vmem=2G -q ${q}
        
        if [[ ${genes_bed} != "" ]]; then
        #intersect the provided genes list with our stats
        echo "intersect the provided genes list with our stats..."
        echo "bedtools intersect -a ${genes_bed} -b ${het_out}/${chr}.ALL.HET.bed -wo > ${gene_stats_folder}/${chr}.genes_HET.bed" |qsub -o ${het_logs}/chr${chr}_\$JOB_ID_bed_intersect.log -e ${het_logs}/chr${chr}_\$JOB_ID_bed_intersect.e -V -N het_bed_intersect_chr${chr} -hold_jid het_concat_chr${chr} -l h_vmem=2G -q ${q}
        fi        

    ;;
    R)

        # STEP 4: ROH stats
        ROH_out="${out_folder}/ROH"
        ROH_logs="${out_folder}/LOGS/ROH"

        mkdir -p ${ROH_out}
        mkdir -p ${ROH_logs}

        if [[ ${genes_bed} != "" ]]; then
        #intersect the provided genes list with our stats
        echo "intersect the provided genes list with our stats..."
        bedtools intersect -a ${genes_bed} -b ${ROH_out}/${chr}.ALL.ROH.bed -wo > ${gene_stats_folder}/${chr}.genes_ROH.bed
        fi

    ;;
    S)
        # echo "/home/cocca/scripts/bash_scripts/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}"
        ${base_bash_scripts}/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}
        # STEP 5: Singleton density
        singletons_out="${out_folder}/singletons"
        singletons_logs="${out_folder}/LOGS/singletons"

        mkdir -p ${singletons_out}
        mkdir -p ${singletons_logs}

        # We need to pass the samples list to the singletons script here
        if [[ ${sample_list} != "" ]]; then
            s_num=$(wc -l ${sample_list} | cut -f 1 -d " ")
            echo " You provided the following sample list: ${sample_list}, containing ${s_num} samples. "
            a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_SINGLETONS.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${singletons_out} ${sample_list}"|qsub -t 1-${a_size} -o ${singletons_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${singletons_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N singletons_chr${chr} -l h_vmem=${m} -q ${q}
        else
            echo " No samples list provided. "
            a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_SINGLETONS.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${singletons_out}"|qsub -t 1-${a_size} -o ${singletons_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${singletons_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N singletons_chr${chr} -l h_vmem=${m} -q ${q}
        fi

        # stats concat step
        echo "${base_bash_scripts}/VPP_GR2018_pipeline_CONCAT.sh -S ${singletons_out} ${chr} ${singletons_out}"|qsub -o ${singletons_logs}/chr${chr}_\$JOB_ID_concat.log -e ${singletons_logs}/chr${chr}_\$JOB_ID_concat.e -V -N singletons_concat_chr${chr} -hold_jid singletons_chr${chr} -l h_vmem=2G -q ${q}

        if [[ ${genes_bed} != "" ]]; then
        #intersect the provided genes list with our stats
        echo "intersect the provided genes list with our stats..."
        echo "bedtools intersect -a ${genes_bed} -b ${singletons_out}/${chr}.ALL.singletons.bed -wo > ${gene_stats_folder}/${chr}.genes_singletons.bed" | qsub -o ${singletons_logs}/chr${chr}_\$JOB_ID_bed_intersect.log -e ${singletons_logs}/chr${chr}_\$JOB_ID_bed_intersect.e -V -N singletons_bed_intersect_chr${chr} -hold_jid singletons_concat_chr${chr} -l h_vmem=2G -q ${q}
        fi
    ;;
    I)

        # STEP 6: IHS
        IHS_out="${out_folder}/IHS"
        IHS_logs="${out_folder}/LOGS/IHS"

        mkdir -p ${IHS_out}
        mkdir -p ${IHS_logs}

        # a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_IHS.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${IHS_out}"|qsub -t 1-${a_size} -o ${IHS_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${IHS_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N IHS_chr${chr} -l h_vmem=${m} -q ${q}
        # echo "${base_bash_scripts}/VPP_GR2018_pipeline_IHS.sh ${chr} ${pop_vcf} ${IHS_out} ${out_folder}"|qsub -o ${IHS_logs}/chr${chr}_\$JOB_ID.log -e ${IHS_logs}/chr${chr}_\$JOB_ID.e -V -N IHS_chr${chr} -l h_vmem=${m} -q ${q}
        # echo "${base_bash_scripts}/VPP_GR2018_pipeline_IHS.sh ${chr} ${pop_vcf} ${IHS_out} ${out_folder} ${m} ${q}"|qsub -o ${IHS_logs}/chr${chr}_\$JOB_ID.log -e ${IHS_logs}/chr${chr}_\$JOB_ID.e -V -N IHS_chr${chr} -l h_vmem=2G -q ${q}
        ${base_bash_scripts}/VPP_GR2018_pipeline_IHS.sh ${chr} ${pop_vcf} ${IHS_out} ${out_folder} ${m} ${q}
        
        # # here we need now to submit a job that will split by 50Kb windows the results
        # # we will use the chunk file generated at the beginning to extract the snps belonging to each windows and the relative iHS score
        # # than we will calculate the average iHS by window and put all back together
        # a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_IHS_chunkinator.sh ${out_folder}/${chr}_chunk_file.txt ${IHS_out}/${chr}.vcf ${IHS_out}"|qsub -t 1-${a_size} -o ${IHS_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${IHS_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N IHS_split_chr${chr} -l h_vmem=${m} -q ${q} -hold_jid IHS_chr${chr} IHS1_chr${chr} IHS2_chr${chr} IHS3_chr${chr} 

        # stats concat step
        echo "${base_bash_scripts}/VPP_GR2018_pipeline_CONCAT.sh -I ${IHS_out} ${chr} ${IHS_out}"|qsub -o ${IHS_logs}/chr${chr}_\$JOB_ID_concat.log -e ${IHS_logs}/chr${chr}_\$JOB_ID_concat.e -V -N IHS_concat_chr${chr} -hold_jid IHS*_chr${chr} -l h_vmem=2G -q ${q}
        # echo "${base_bash_scripts}/VPP_GR2018_pipeline_CONCAT.sh -I ${IHS_out} ${chr} ${IHS_out}"|qsub -o ${IHS_logs}/chr${chr}_\$JOB_ID_concat.log -e ${IHS_logs}/chr${chr}_\$JOB_ID_concat.e -V -N IHS_concat_chr${chr} -l h_vmem=2G -q ${q}
        if [[ ${genes_bed} != "" ]]; then
        #intersect the provided genes list with our stats
        echo "intersect the provided genes list with our stats..."
        echo "bedtools intersect -a ${genes_bed} -b ${IHS_out}/${chr}.ALL.IHS.bed -wo > ${gene_stats_folder}/${chr}.genes_IHS.bed" | qsub -o ${IHS_logs}/chr${chr}_\$JOB_ID_bed_intersect.log -e ${IHS_logs}/chr${chr}_\$JOB_ID_bed_intersect.e -V -N IHS_bed_intersect_chr${chr} -hold_jid IHS_concat_chr${chr} -l h_vmem=2G -q ${q}
        fi

    ;;
    C)
        # echo "/home/cocca/scripts/bash_scripts/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}"
        ${base_bash_scripts}/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}
        # STEP 7: CADD score distribution
        CADD_out="${out_folder}/CADD"
        CADD_logs="${out_folder}/LOGS/CADD"

        mkdir -p ${CADD_out}
        mkdir -p ${CADD_logs}

        a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_CADD.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${CADD_out}"|qsub -t 1-${a_size} -o ${CADD_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${CADD_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N CADD_chr${chr} -l h_vmem=${m} -q ${q}

        # stats concat step
        echo "${base_bash_scripts}/VPP_GR2018_pipeline_CONCAT.sh -C ${CADD_out} ${chr} ${CADD_out}"|qsub -o ${CADD_logs}/chr${chr}_\$JOB_ID_concat.log -e ${CADD_logs}/chr${chr}_\$JOB_ID_concat.e -V -N CADD_concat_chr${chr} -hold_jid CADD_chr${chr} -l h_vmem=2G -q ${q}

        if [[ ${genes_bed} != "" ]]; then
        #intersect the provided genes list with our stats
        echo "intersect the provided genes list with our stats..."
        echo "bedtools intersect -a ${genes_bed} -b ${CADD_out}/${chr}.ALL.CADD.bed -wo > ${gene_stats_folder}/${chr}.genes_CADD.bed" | qsub -o ${CADD_logs}/chr${chr}_\$JOB_ID_bed_intersect.log -e ${CADD_logs}/chr${chr}_\$JOB_ID_bed_intersect.e -V -N CADD_bed_intersect_chr${chr} -hold_jid CADD_concat_chr${chr} -l h_vmem=2G -q ${q}
        fi

    ;;
esac
done
