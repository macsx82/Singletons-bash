#!/usr/bin/env bash

#script to calculate singleton density based in coding and non coding regions

ref_seq=$1
cds_uniq=$2
# sing_file=/shared/Singleton_Boost_PJ/singleton_score/VBI/VBI_${chr}_ALL.singletons
sing_file=$3
outfile=$4
s_mode=$5 #alternatives to calculate singleton scores: POP to use standard population mode, SAMPLE to use custom single sample mode to calculate singleton stats by sample (to be used for association purposes)

#we need to use the Gene ID because it's uniq, and we can use it to link CDS back to the correct gene, even if we have some duplicates in GENCODE data, due to HAVANA/ENSEMBL discrepancies
#if there is  no bed file for the singleton info, we will generate it
if [[ ! -s ${sing_file}.bed ]]; then
echo "Generating bed files for singletons..."
# fgrep -v "INDV" ${sing_file} |awk '{OFS="\t"}{print "chr"$1,$2,$2+length($4),$3,$5}' > ${sing_file}.bed
fgrep -v "INDV" ${sing_file} |awk '{OFS="\t"}{print "chr"$1,$2-1,$2+length($4)-1,$3,$5}' > ${sing_file}.bed
fi

#activate python env needed
source activate py36

case ${s_mode} in
	POP )
		${MY_BASH_SCRIPTS}/singletons_score.py --sing ${sing_file}.bed --gen_seq ${ref_seq} --cds_seq ${cds_uniq} --out ${outfile}

		#Add final singleton scores calculation
		Rscript --no-save --verbose ${MY_BASH_SCRIPTS}/Function_singleton_calculation.r ${outfile} ${outfile}.scores

	;;
	SAMPLE )
		# here we need to extract data for each sample, so we need to cicle throug the singleton bed file and generate a singleton file for each sample
		# Than we can run the score calculation code for each sample and generate a file for each sample and than a summary collecting all data for all samples
		# This is intended to be run by gene
		genome_legth=$6 #Genome length for genome wide singleton density calculation
		# genome_legth=3095677412
		samples=$(cut -f 5 ${sing_file}.bed | sort | uniq )

		#generate a bed file for each sample and run the singleton score calculation
		for sample in ${samples}
		do
			outpath=$(dirname ${outfile})
			s_sing_file_name=$(basename ${sing_file})
			fgrep -w ${sample} ${sing_file}.bed > ${outpath}/${s_sing_file_name}_${sample}.bed

			#get GW singleton number fot the current sample
			s_sing_num=$(wc -l ${outpath}/${s_sing_file_name}_${sample}.bed| cut -f 1 -d " ")


			${MY_BASH_SCRIPTS}/singletons_score.py --sample --sample_id ${sample} --genome_length ${genome_legth} --s_sing_num ${s_sing_num} --sing ${outpath}/${s_sing_file_name}_${sample}.bed --gen_seq ${ref_seq} --cds_seq ${cds_uniq} --out ${outfile}_${sample}.txt
			#before cleanup, I want the singleton density for that sample, but genome wide
			#So we need to extract the singleton number from ALL chromosomes

			# clean up after score generation
			rm ${outpath}/${s_sing_file_name}_${sample}.bed
		done
		
	;;
	*)
		echo "Accepted options for s_mode parameter are POP (to calculate population wide stats) or SAMPLE (to calculate sample levels stats"
	;;
esac

echo "Done singleton score calculation"