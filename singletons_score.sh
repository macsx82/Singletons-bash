#!/usr/bin/env bash

#script to calculate singleton density based in coding and non coding regions
base_script_folder=$(dirname $0)
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
# fgrep -v "INDV" ${sing_file} |awk '{OFS="\t"}{print "chr"$1,$2-1,$2+length($4)-1,$3,$5}' > ${sing_file}.bed
fgrep -v "INDV" ${sing_file} |awk '{OFS="\t"}{print $1,$2-1,$2+length($4)-1,$3,$5}' > ${sing_file}.bed
fi

#activate python env needed
source activate py36

case ${s_mode} in
	POP )
		${base_script_folder}/singletons_score.py --sing ${sing_file}.bed --gen_seq ${ref_seq} --cds_seq ${cds_uniq} --out ${outfile}

		#Add final singleton scores calculation
		Rscript --no-save --verbose ${base_script_folder}/Function_singleton_calculation.r ${outfile} ${outfile}.scores

	;;
	SAMPLE )
		#WORK IN PROGRESS
	;;
	*)
		echo "Accepted options for s_mode parameter are POP (to calculate population wide stats) or SAMPLE (to calculate sample levels stats"
	;;
esac

echo "Done singleton score calculation"