#!/usr/bin/env bash


###############################################
#23/04/2018
#
#Wrapper script for singletons scores calculation
# v0.1
###############################################

# Pipeline workflow:

# Singleton density


#TODO:define a function to generate bed files with region boundaries
#need to find the updated coordinates for each chr, since the absolute dist
base_bash_scripts=$(dirname $0)
######################
# ARGS:
# $1 is reserved for pipeline options
#chr number (TODO: parse chr name in vcf to check if the chr name is just the number or chrN)
chr=$2
#size in bp
win_size=$3 
out_folder=$4
pop_vcf=$5
genome_build=$6
m=$7
q=$8
#we need to account for a sample list,and for a bed file containing genes
#The list will be an INCLUSION list, with a single column of samples's ids to include in the analyses
sample_list=$9
#this bed file is useless, most of the times, we need to check if it's still usefull to have this option
genes_bed=${10}
# pipe_step=$7
# STEP 1: create a bed file for each chromosome to work with regions in parallel
# a) retrieve updated chromosome coordinates
case ${genome_build} in
    GRCh37 )
        chrom=(`mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B --skip-column-names -e "select chrom, 0, size as coords  from hg19.chromInfo where chrom = 'chr${chr}';"`)
        ;;
    GRCh38 )
        chrom=(`mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B --skip-column-names -e "select chrom, 0, size as coords  from hg38.chromInfo where chrom = 'chr${chr}';"`)
        ;;
esac
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
while getopts ":THRSICh" opt; do
    case ${opt} in
    S)
        ${base_bash_scripts}/generic_chunk_generator.sh ${chunk_mode} ${start}-${end} ${win_size} ${chr} ${chunk_outfile}
        # Singleton density
        singletons_out="${out_folder}/singletons"
        singletons_logs="${out_folder}/LOGS/singletons"

        mkdir -p ${singletons_out}
        mkdir -p ${singletons_logs}

        # We need to pass the samples list to the singletons script here
        if [[ ${sample_list} != "" ]]; then
            s_num=$(wc -l ${sample_list} | cut -f 1 -d " ")
            echo " You provided the following sample list: ${sample_list}, containing ${s_num} samples. "
            #generating a job array to work by chunk: this is hard coded for SGE, but more options are coming
            a_size=`wc -l ${out_folder}/${chr}_chunk_file.txt| cut -f 1 -d " "`;echo "${base_bash_scripts}/ja_runner_par_TRST.sh -t ${base_bash_scripts}/VPP_GR2018_pipeline_SINGLETONS.sh ${out_folder}/${chr}_chunk_file.txt ${pop_vcf} ${singletons_out} ${sample_list}"|qsub -t 1-${a_size} -o ${singletons_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${singletons_logs}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N singletons_chr${chr} -l h_vmem=${m} -q ${q}
        else
            echo " No samples list provided. "
            #generating a job array to work by chunk: this is hard coded for SGE, but more options are coming
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
    h)
    echo "Script usage:"
    echo "VPP_GR2018_pipeline.sh [workflow option] [arguments] "
    echo "workflow options:"
    echo "                 -S: singletons workflow "
    echo ""
    echo "arguments: "
    echo "In this preliminary version of the code, all arguments are positional:"
    echo -e "$2 : chromosome \n
         $3 : window size (in bp)  \n
         $4 : output folder \n
         $5 : vcf file \n
         $6 : genome_build \n
         $7 : memory requirement for job submission \n
         $8 : queue for job submission \n
         $9 : sample INCLUSION list, with a single column of samples's ids to include in the analyses \n
         $10 : bed file containing genes [mostly deprecated] \n
        "
    ;;
    
esac
done
