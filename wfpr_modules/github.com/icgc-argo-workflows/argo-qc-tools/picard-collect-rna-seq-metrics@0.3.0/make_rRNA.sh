#!/usr/bin/env bash
# make_rRNA.sh
# Kamil Slowikowski
# December 12, 2014
#
# Modified: Arindam Ghosh (July 24, 2019)
#           Linda Xiang (August 25, 2021)
#
#
# Referenc Genome: GRCh38
# 1. Prepare chromosome sizes file from fasta sequence if needed.
#     E.g., download fasta sequence from ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#     samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
#     cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > sizes.genome
#
# 2. Make an interval_list file suitable for CollectRnaSeqMetrics.jar using Gencode gene annotation
# E.g., http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
#
# Picard Tools CollectRnaSeqMetrics.jar:
#   https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics

# Chromosome sizes file
chrom_sizes=$1

# Genes annotation from Gencode.
genes=$2

# Output file suitable for Picard CollectRnaSeqMetrics.jar.
rRNA=GRCh38.rRNA.interval_list

# Sequence names and lengths. (Must be tab-delimited.)
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' $chrom_sizes | \
    grep -v _ \
> $rRNA

# Intervals for rRNA transcripts.
grep 'gene_type "rRNA_pseudogene"' $genes | \
    awk '$3 == "gene"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /gene_id "([^"]+)"/ or die "no gene_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> $rRNA