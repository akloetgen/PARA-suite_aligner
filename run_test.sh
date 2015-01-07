#!/bin/bash


echo "run bwa alnep"
echo "bwa alnep -p test/error_profile_MSI1+MSI2_bowtie2.tsv -t 12 -f test/testout_ep.sai /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa /home/akloetgen/read_mapper/benchmark/benchmark_refseq_ensembl75_own-100T2C.fastq"

bwa alnep -p test/error_profile_MSI1+MSI2_bowtie2.tsv -t 2 -f test/testout_ep.sai /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa test/1read_mut.fastq
#bwa alnep -p test/error_profile_MSI1+MSI2_bowtie2.tsv -n 2 -t 20 -f test/testout_ep.sai /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa /home/akloetgen/read_mapper/benchmark/benchmark_refseq_ensembl75_own-100T2C.fastq
#bwa alnep -p test/error_profile_MSI1+MSI2_bowtie2.tsv -n 2 -t 20 -f test/testout_ep2.sai /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa /net/cancergenomics/projects/parclip_andreas/data/PARCLIP_test/MSI2/149_MSI2_lane6.trimmed4.fastq
#bwa aln -n 2 -t 20 -f test/testout.sai /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa /home/akloetgen/read_mapper/benchmark/benchmark_refseq_ensembl75_own-100T2C.fastq

bwa samse -f test/testout_ep.sam /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa test/testout_ep.sai test/1read_mut.fastq
#bwa samse -f test/testout_ep.sam /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa test/testout_ep.sai /home/akloetgen/read_mapper/benchmark/benchmark_refseq_ensembl75_own-100T2C.fastq
#bwa samse -f test/testout_ep2.sam /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa test/testout_ep2.sai /net/cancergenomics/projects/parclip_andreas/data/PARCLIP_test/MSI2/149_MSI2_lane6.trimmed4.fastq
#bwa samse -f test/testout.sam /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa test/testout.sai /home/akloetgen/read_mapper/benchmark/benchmark_refseq_ensembl75_own-100T2C.fastq

echo "bam conv"
samtools view -bS -t /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.8/hg19.fa test/testout_ep.sam -o test/testout_ep.bam
echo "filter"
samtools view -q 10 -b test/testout_ep.bam -o test/testout_ep.unique.bam
echo "sort"
samtools sort test/testout_ep.unique.bam test/testout_ep.unique.sort
echo "index"
samtools index test/testout_ep.unique.sort.bam
