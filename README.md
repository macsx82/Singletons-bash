# Singletons-bash
Bash scripts needed for singleton scores and stats calculation

This bash version of the pipeline is based on the work presented in:

Mezzavilla M, Cocca M, Guidolin F, Gasparini P. A population-based approach for gene prioritization in understanding complex traits. Hum Genet. 2020;139(5):647-655. doi:10.1007/s00439-020-02152-4

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
	* This step will be performed in any case, so it should be a stand-alone procedure, to give the option to reuse already generated data

3. Calculate singleton stats
	+ This step will allow two option:
		1. Genomic region stats
			+ Singleton density for each gene across all samples
			+ Singleton count for each gene across all samples
		2. Sample level stats
			+ Genome wide singleton density for each sample
			+ Singleton density for each sample in the provided genomic region
			+ Singleton count for each sample in the provided genomic region

4. Calculate singleton scores
	+ This step will be performed only in case of genomic region pipeline

---
##Input files

+ Multisample vcf files for the study population
+ ...

---
##Sample command

Generate data files for b38 on ORFEO

```bash
zcat /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation.gff3.gz| grep -E -v "^#"| awk '$3=="gene"'| fgrep "protein_coding" > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation.gff3.GENE.bed
zcat /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation.gff3.gz| grep -E -v "^#"| awk '$3=="CDS"'| fgrep "transcript_type=protein_coding" > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation.gff3.CDS.bed

cut -f 1,4,5,9 /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation.gff3.GENE.bed | awk 'BEGIN{OFS="\t"}{split($4,a,";");split(a[2],g,"=");split(a[4],n,"="); print $1,$2,$3,n[2],g[2]}' > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation_GENES.bed
cut -f 1,4,5,9 /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation.gff3.CDS.bed | awk 'BEGIN{OFS="\t"}{split($4,a,";");split(a[3],g,"=");split(a[6],n,"=");split(a[2],p,"=");split(a[4],t,"="); print $1,$2,$3,n[2],g[2],p[2],t[2]}' > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v34.annotation_CDS.bed

zcat /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation.gff3.gz| grep -E -v "^#"| awk '$3=="gene"'| fgrep "protein_coding" > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation.gff3.GENE.bed
zcat /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation.gff3.gz| grep -E -v "^#"| awk '$3=="CDS"'| fgrep "transcript_type=protein_coding" > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation.gff3.CDS.bed

cut -f 1,4,5,9 /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation.gff3.GENE.bed | awk 'BEGIN{OFS="\t"}{split($4,a,";");split(a[2],g,"=");split(a[4],n,"="); print $1,$2,$3,n[2],g[2]}' > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation_GENES.bed
cut -f 1,4,5,9 /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation.gff3.CDS.bed | awk 'BEGIN{OFS="\t"}{split($4,a,";");split(a[3],g,"=");split(a[6],n,"=");split(a[2],p,"=");split(a[4],t,"="); print $1,$2,$3,n[2],g[2],p[2],t[2]}' > /storage/burlo/cocca/resources/hgRef/hg38/gencode.v35.annotation_CDS.bed


for gencode_v in 34 35
do 
ref_seq=/storage/burlo/cocca/resources/hgRef/hg38/gencode.v${gencode_v}.annotation_GENES.bed
cds_uniq=/storage/burlo/cocca/resources/hgRef/hg38/gencode.v${gencode_v}.annotation_CDS.bed
outfolder=/storage/burlo/cocca/resources/hgRef/hg38
mkdir -p /storage/burlo/cocca/resources/hgRef/hg38/CHR
for chr in {1..22} X Y M
do
# chr=1
    awk -v c=${chr} '$1=="chr"c' ${ref_seq} > ${outfolder}/gencode.v${gencode_v}.${chr}.annotation_GENES.bed
    cut -f 1,5 ${outfolder}/gencode.v${gencode_v}.${chr}.annotation_GENES.bed > ${outfolder}/gencode.v${gencode_v}.${chr}.GENES.bed
    awk -v c=${chr} '$1=="chr"c' ${cds_uniq} > ${outfolder}/gencode.v${gencode_v}.${chr}.annotation_CDS.bed
done
done
```

We need to run the singleton scores calculation for the INGI data lifted to the GRCh38 build

```bash
for pop in CAR FVG VBI
do

for chr in {1..22} X
do

#chr=22
win_size=50000
out_folder=/home/cocca/analyses/SINGLETONS/20201021/${pop}/${chr}
pop_vcf=/shared/INGI_WGS/GRCh38/${pop}/sorted/${chr}.ucsc_sorted.vcf.gz
m=5G
q="fast"

/home/cocca/scripts/bash_scripts/VPP_GR2018_pipeline.sh -S ${chr} ${win_size} ${out_folder} ${pop_vcf} ${m} ${q}

done
done

```

Calculate singleton scores

```bash
for gencode_v in 34 35
do

for pop in CAR FVG VBI
do
for chr in {1..22} X   
do
#list of uniq gene ids, by chr
ref_seq=/shared/resources/hgRef/hg38/CHR/gencode.v${gencode_v}.${chr}.annotation_GENES.bed
cds_uniq=/shared/resources/hgRef/hg38/CHR/gencode.v${gencode_v}.${chr}.annotation_CDS.bed

sing_file=/home/cocca/analyses/SINGLETONS/20201021/${pop}/${chr}/singletons/${chr}.ALL.samples.singletons
out_path=/home/cocca/analyses/SINGLETONS/20201021/${pop}/gencode_v${gencode_v}/singleton_score
outfile=${out_path}/${pop}_SING_SCORE.${chr}.txt
log_folder=${out_path}/logs

job_name=s_score_${pop}_${chr}
qu="fast,all.q"

mkdir -p ${log_folder}

echo "~/scripts/bash_scripts/singletons_score.sh ${ref_seq} ${cds_uniq} ${sing_file} ${outfile} POP"|qsub -o ${log_folder}/\$JOB_ID_${pop}_${chr}.log -e ${log_folder}/\$JOB_ID_${pop}_${chr}.e -V -N ${job_name} -l h_vmem=4G -q ${qu}

done
done
done
