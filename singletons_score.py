#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Massimiliano Cocca
:contact: massimiliano.cocca [at] burlo.trieste.it
:Date: *20.03.2019

:Description:

Calculate singleton score by gene and CDS based on gene_id and transcript_id

"""

__author__ = 'macsx82'

import argparse
import gzip
import re
import sys
import collections
import numpy
import subprocess


# args.infile="/shared/Singleton_Boost_PJ/singleton_score/VBI/VBI_1_ALL.singletons.bed"
# args.gen_seq="/shared/Singleton_Boost_PJ/ref_data/gencode.v29lift37.1.annotation_GENES.bed"
# args.cds_seq="/shared/Singleton_Boost_PJ/ref_data/gencode.v29lift37.1.annotation_CDS.bed"
# args.outfile="/shared/Singleton_Boost_PJ/singleton_score/VBI/SING_SCORE.1.txt"


def main(args):

    # sing_f = gzip.open(args.infile) if args.infile.endswith('.gz') else open(args.infile)
    gene_f = gzip.open(args.gen_seq) if args.gen_seq.endswith('.gz') else open(args.gen_seq)
    cds_f = gzip.open(args.cds_seq) if args.cds_seq.endswith('.gz') else open(args.cds_seq)
    # gnomad_f = gzip.open(args.gnomad_seq) if args.gnomad_seq.endswith('.bgz') else open(args.gnomad_seq)
    # o_tab = gzip.open(args.outfile, 'w') if args.outfile.endswith('.gz') else open(args.outfile,'w')

    #generate the intersection using bedtools as in the shell script for genes and cds
    cmd_gene = 'bedtools intersect -a ' + args.infile +' -b '+ args.gen_seq + ' -wo | cut -f 10 | sort | uniq -c| awk \'{OFS="\\t"}{print $1,$2,$3}\''
    cmd_cds = 'bedtools intersect -a ' + args.infile +' -b '+ args.cds_seq + ' -wo | cut -f 10,12 | sort | uniq -c| awk \'{OFS="\\t"}{print $1,$2,$3}\''

    # n_s_tot=subprocess.run(cmd_1,stdout=subprocess.PIPE,shell=True)
    n_s_tot=subprocess.run(cmd_gene,stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').split("\n")
    cds_s_tot=subprocess.run(cmd_cds,stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').split("\n")

    #need to transform all the singletons count into in dictionaries to combine with the others
    n_s_tot_list=collections.defaultdict(list)
    for s_count_g in n_s_tot:
        if s_count_g != "":
            s_count_g=s_count_g.strip().split("\t")
            gene_id=s_count_g[1].split(".")[0]
            s_count=s_count_g[0]
            n_s_tot_list[gene_id].append(int(s_count))
            # n_s_tot_list[gene_id]=[gene_name], int(s_count)]

    cds_s_tot_list=collections.defaultdict(lambda: collections.defaultdict(list))
    for s_count_cds in cds_s_tot:
        if s_count_cds != "":
            s_count_cds=s_count_cds.strip().split("\t")
            gene_id=s_count_cds[1].split(".")[0]
            transcript_id=s_count_cds[2].split(".")[0]
            s_count=s_count_cds[0]
            cds_s_tot_list[gene_id][transcript_id].append(int(s_count))

    #calculate gene length and create a dictionary
    gene_list={gene.strip().split("\t")[4].split(".")[0]:[gene.strip().split("\t")[3],int(gene.strip().split("\t")[2])-int(gene.strip().split("\t")[1])+1] for gene in gene_f }
    # gene_list={}
    # for gene in gene_f:
    #     gene=gene.strip().split("\t")
    #     gene_id=gene[4]
    #     gene_length=int(gene[2])-int(gene[1])+1
    #     gene_list[gene_id]=gene_length

    #now we need the same for each CDS, but we need a dictionary with gene_id AND transcript id, for the size of each cds region and for the total cds size for each transcript
    cds_list=collections.defaultdict(lambda: collections.defaultdict(list))
    for cds in cds_f:
        cds=cds.strip().split("\t")
        gene_id=cds[4].split(".")[0]
        transcript_id=cds[6].split(".")[0]
        cds_length=int(cds[2])-int(cds[1])+1
        cds_list[gene_id][transcript_id].append(cds_length)

    
    if args.sample:
        #We are working in sample mode, so we want to calculate wg density
        #Calculate Genome wide singleton density for a single sample mode
        gw_s_density=float(args.s_sing_num)/float(args.genome_length)

        # In this case, we want to generate two separate files 
        # gnomad_f = gzip.open(args.gnomad_seq) if args.gnomad_seq.endswith('.bgz') else open(args.gnomad_seq)
        gene_f_out=args.outfile + "_GENE"
        cds_f_out=args.outfile + "_CDS"
        gene_out = gzip.open(gene_f_out, 'w') if gene_f_out.endswith('.gz') else open(gene_f_out,'w')
        cds_out = gzip.open(cds_f_out, 'w') if cds_f_out.endswith('.gz') else open(cds_f_out,'w')

        #Last step is to generate the report going through all genes and transcripts and CDS
        #add header
        header_line_GENE=["SAMPLE_ID","GENE_ID","GENE_NAME","SING_GENE","GENE_LENGTH","GENE_SING_DENSITY","GW_SINGLETON_DENSITY","NORMALIZED_SINGLETON_DENSITY"]
        header_line_CDS=["SAMPLE_ID","GENE_ID","GENE_NAME","TRANSCRIPT_ID","SING_TRANSCRIPT","TRANSCRIPT_LENGTH","CDS_SING_DENSITY","GW_SINGLETON_DENSITY","NORMALIZED_SINGLETON_DENSITY"]
        
        print('%s' %("\t".join(header_line_GENE)), file=gene_out)
        print('%s' %("\t".join(header_line_CDS)), file=cds_out)
        
        for gene,g_length in gene_list.items():
            #get the singleton count by gene
            try:
                count_s_gene=sum(n_s_tot_list[gene])
            except:
                #if there are no overlap with any gene, return 0
                count_s_gene=0
            g_s_density=float(count_s_gene)/float(g_length[1])
            #now print the report by gene
            splitted_line_gene=[args.sample_name,gene,g_length[0],count_s_gene,g_length[1],g_s_density,gw_s_density,float(g_s_density)/float(gw_s_density)]
            print('%s' %("\t".join(list(map(str,splitted_line_gene)))), file=gene_out)
            
            #now get the count for each transcript of this gene
            for transcript,t_length in cds_list[gene].items():
                try:
                    curr_transcript_s_count=sum(cds_s_tot_list[gene][transcript])
                except:
                    #if there are no overlap with any transcript, return 0
                    curr_transcript_s_count=0
                cds_s_density=float(curr_transcript_s_count)/float(sum(t_length))
                #now print the report by transcript
                splitted_line_cds=[args.sample_name,gene,g_length[0],transcript,curr_transcript_s_count,sum(t_length),cds_s_density,gw_s_density,float(cds_s_density)/float(gw_s_density)]
                # print(splitted_line)
                print('%s' %("\t".join(list(map(str,splitted_line_cds)))), file=cds_out)

        gene_f.close()
        cds_f.close()
        gene_out.close()
        cds_out.close()

    else :

        # gnomad_f = gzip.open(args.gnomad_seq) if args.gnomad_seq.endswith('.bgz') else open(args.gnomad_seq)
        o_tab = gzip.open(args.outfile, 'w') if args.outfile.endswith('.gz') else open(args.outfile,'w')

        #Last step is to generate the report going through all genes and transcripts and CDS
        #add header
        header_line=["GENE_ID","GENE_NAME","SING_GENE","GENE_LENGTH","GENE_SING_DENSITY","TRANSCRIPT_ID","SING_TRANSCRIPT","TRANSCRIPT_LENGTH","CDS_SING_DENSITY"]
        print('%s' %("\t".join(header_line)), file=o_tab)

        for gene,g_length in gene_list.items():
            #get the singleton count by gene
            try:
                count_s_gene=sum(n_s_tot_list[gene])
            except:
                #if there are no overlap with any gene, return 0
                count_s_gene=0
            g_s_density=float(count_s_gene)/float(g_length[1])
            #now get the count for each transcript of this gene
            for transcript,t_length in cds_list[gene].items():
                try:
                    curr_transcript_s_count=sum(cds_s_tot_list[gene][transcript])
                except:
                    #if there are no overlap with any transcript, return 0
                    curr_transcript_s_count=0
                cds_s_density=float(curr_transcript_s_count)/float(sum(t_length))
                #now print the report
                splitted_line=[gene,g_length[0],count_s_gene,g_length[1],g_s_density,transcript,curr_transcript_s_count,sum(t_length),cds_s_density]
                # print(splitted_line)
                print('%s' %("\t".join(list(map(str,splitted_line)))), file=o_tab)
                # print >> o_tab, '%s' %("\t".join(list(map(str,splitted_line))))

        gene_f.close()
        cds_f.close()
        o_tab.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #define all input parameteres and files
    parser.add_argument('--sing', dest="infile", help='Input singleton file in bed format')
    parser.add_argument('--gen_seq', dest="gen_seq", help='Gene regions from GENCODE in bed format')
    parser.add_argument('--cds_seq', dest="cds_seq", help='CDS regions from GENCODE in bed format')
    parser.add_argument('--gnomad', dest="gnomad_seq", help='GnomAD data by transcript in bgz format')
    parser.add_argument('--sample', action="store_true", help='Flag to select the single sample mode')
    parser.add_argument('--sample_id', dest="sample_name", help='Sample id under analysis')
    parser.add_argument('--s_sing_num', dest="s_sing_num", help='Singleton number per sample')
    parser.add_argument('--genome_length', dest="genome_length", help='Provided Genome length to calculate genome wide singleton density by sample')
    parser.add_argument("-o","--out", dest="outfile", help="Output file path",default="sing_score.tab")
    args = parser.parse_args()

print(args)
main(args)

