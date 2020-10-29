# Singletons-bash
Bash scripts needed for singleton scores and stats calculation

This bash version of the pipeline is based on the work presented in:

Mezzavilla M, Cocca M, Guidolin F, Gasparini P. A population-based approach for gene prioritization in understanding complex traits. Hum Genet. 2020;139(5):647-655. doi:10.1007/s00439-020-02152-4

---
#DISCLAIMER

This is the very first version of the pipeline, meant to be shared quickly while we work on a more polished and portable Snakemake version.

There are some hard coded parameters you should modify if you plan to run the scripts on a cluster not based on the SGE queue system.
	
* in the folder there is a ja_runner_par.sh script to be used on LSB cluster: you can use this script to modify the job array submission syntax in the VPP_GR2018_pipeline.sh script (lines 78 and 82)
* you should modify job submission syntax according to your cluster on lines 86 and 91 in VPP_GR2018_pipeline.sh


---
##Requirements

+ Python <= 3.6
+ bedtools
+ vcftools
+ bcftools

---

##Pipeline steps:

1. Provide genomic region input
	* Get gene info from GENCODE data (if requested)
	* Provide a genomic region using chr:start-end informations

2. Calculate singletons for the study population
	* Computationally intensive step
	* This step will be performed in any case, and once per population, so the data can be reused

3. Calculate singleton stats
	+ This step will allow two option:
		1. Genomic region stats
			+ Singleton density for each gene across all samples
			+ Singleton count for each gene across all samples
			+ Singleton scores (SSC and DSC)
			
		2. Sample level stats (in development )
			+ Genome wide singleton density for each sample
			+ Singleton density for each sample in the provided genomic region
			+ Singleton count for each sample in the provided genomic region

---
##Input files

+ Multisample vcf files for the study population
+ Bed files for genomic regions comprising whole gene regions and coding regions only (examples provided to extract the correct data format from GENCODE data)

---
##Sample command

Generate genomic regions bed files for b38 using GENCODE v35 data:

* Extract only genes annotated as protein coding:

```bash
zcat [..]/hg38/gencode.v35.annotation.gff3.gz| grep -E -v "^#"| awk '$3=="gene"'| fgrep "protein_coding" > [..]/hg38/gencode.v35.annotation.gff3.GENE.bed
zcat [..]/hg38/gencode.v35.annotation.gff3.gz| grep -E -v "^#"| awk '$3=="CDS"'| fgrep "transcript_type=protein_coding" > [..]/hg38/gencode.v35.annotation.gff3.CDS.bed
```

* Generate the final formatted bed files:

```bash
cut -f 1,4,5,9 [..]/hg38/gencode.v35.annotation.gff3.GENE.bed | awk 'BEGIN{OFS="\t"}{split($4,a,";");split(a[2],g,"=");split(a[4],n,"="); print $1,$2,$3,n[2],g[2]}' > [..]/hg38/gencode.v35.annotation_GENES.bed
cut -f 1,4,5,9 [..]/hg38/gencode.v35.annotation.gff3.CDS.bed | awk 'BEGIN{OFS="\t"}{split($4,a,";");split(a[3],g,"=");split(a[6],n,"=");split(a[2],p,"=");split(a[4],t,"="); print $1,$2,$3,n[2],g[2],p[2],t[2]}' > [..]/hg38/gencode.v35.annotation_CDS.bed
```

* Split bed files by chromosome for parallelization:

```bash
ref_seq="hg38/gencode.v35.annotation_GENES.bed"
cds_uniq="hg38/gencode.v35.annotation_CDS.bed"
outfolder="hg38/CHR"
mkdir -p ${outfolder}
for chr in {1..22} X Y M
do
# chr=1
    awk -v c=${chr} '$1=="chr"c' ${ref_seq} > ${outfolder}/gencode.v${gencode_v}.${chr}.annotation_GENES.bed
    cut -f 1,5 ${outfolder}/gencode.v${gencode_v}.${chr}.annotation_GENES.bed > ${outfolder}/gencode.v${gencode_v}.${chr}.GENES.bed
    awk -v c=${chr} '$1=="chr"c' ${cds_uniq} > ${outfolder}/gencode.v${gencode_v}.${chr}.annotation_CDS.bed
done
done
```


* Run the singleton scores calculation for our study population:

```bash

pop="TEST_POP"

for chr in {1..22} X
do
win_size=50000
out_folder=[..]/SINGLETONS/20201021/${pop}/${chr}
pop_vcf=[base-population-folder]/${pop}/${chr}.vcf.gz

#the code is actually written to run on an SGE cluster, but can be modified to run on all sort of clusters
m=5G #memory requiremente for job submission on a queue system
q="fast" #queue to submit the job

/home/cocca/scripts/bash_scripts/VPP_GR2018_pipeline.sh -S ${chr} ${win_size} ${out_folder} ${pop_vcf} ${m} ${q}

done
done

```

* Calculate singleton scores

```bash
pop="TEST_POP"

for chr in {1..22} X   
do
#list of uniq gene ids, by chr that we generated earlier
ref_seq=../hg38/CHR/gencode.v35.${chr}.annotation_GENES.bed
cds_uniq=../hg38/CHR/gencode.v35.${chr}.annotation_CDS.bed

sing_file=[base-population-folder]/${pop}/${chr}/singletons/${chr}.ALL.samples.singletons #this file come from the previous step
out_path=[base-output-path]/SINGLETONS/20201021/${pop}/singleton_score
outfile=${out_path}/${pop}_SING_SCORE.${chr}.txt
log_folder=${out_path}/logs

job_name=s_score_${pop}_${chr}
#the code is actually written to run on an SGE cluster, but can be modified to run on all sort of clusters
qu="fast,all.q"

mkdir -p ${log_folder}

echo "~/scripts/bash_scripts/singletons_score.sh ${ref_seq} ${cds_uniq} ${sing_file} ${outfile} POP"|qsub -o ${log_folder}/\$JOB_ID_${pop}_${chr}.log -e ${log_folder}/\$JOB_ID_${pop}_${chr}.e -V -N ${job_name} -l h_vmem=4G -q ${qu}

done
done
done
```
