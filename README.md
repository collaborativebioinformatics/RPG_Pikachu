# Telomere-to-telomere genome for annotating common alleles & abnormally avoided stop-codons

<img align="right" width="120" height="200" src="https://github.com/collaborativebioinformatics/popchrom/blob/main/others/pikachu2.jpeg">

## Contributors 
-  Elbay Aliyev - Lead 
-  Muhammad Sohail Raza - Writer
-  Shangzhe Zhang, Anastasia Illarionova, Prashant Ranjan, Tiancheng Xu, Aditi Sammi - Tech support 
-  Bryce Kille, ChunHsuan LO - SysAdmin

## Goal 
- To identify rare variants in CHM13 and recorrecting some of them as common variants. By Checking the fasta file directly.
- To screen out in-frame stop codons sites that disagree with Ribo-seq analysis for validting the annotation.
- Propose a more representative reference genome by 1000 Genomes Project 

## Flowchart

by Anastasia Illarionova

<img width="1200" alt="flowchart" src="https://github.com/collaborativebioinformatics/popchrom/blob/main/others/flowchart_version5.png">

## Introduction 

After two decades of refinements, the human reference genome (GRCh38) has become more accurate and mostly complete. However, there are still hundreds of unresolved gaps persist, and no single chromosome has been finished from end to end because of the existence of highly repeated sequences (called transposable elements). Foutunatedly, by using the high-coverage & ultra-long-read technologies, several scientific groups have presented a new human genome assembly that surpasses the continuity of GRCh38, along with a gapless, telomere-to-telomere assembly of a human chromosome, which is called CHM13 reference genome. 

## Input:

- vcf files (using CHM13_v1.0 as reference already)

## Outputs: 

- New reference Fasta file (for all samples in 1GP)
- (Optional) New reference Fasta files for 5 subpopulations
- ClinVar annotation of common alleles

## Methods

### Dependancies

- python3
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [bcftools](https://samtools.github.io/bcftools/)
- [gffread](https://github.com/gpertea/gffread)
- [GATK](https://gatk.broadinstitute.org)

### Sample Workflow (chr22)

#### I. Data Acquisition and Preprocessing:

1. Downloading the raw CHM13 fasta file.

2. Downloading the raw population based VCF (CHM13 fasta based).

3. Downloading the raw Annotation file (gff3).

#### II. Core tasks:

**1.** Identifying stop codon sites in CHM13 fasta file (Shangzhe, Muhamad, Bryce)

  -  Extract the protein sequences from GFF annotation and FASTA file

  ```
  gffread -g chm13.draft_v1.0.fasta chm13.draft_v1.0.gene_annotation.v4.gff3 -y chm13.v1.0.pep.fasta
  ```

  -  Costom script to extract the positions of in-frame stop codons (Only for chr22)

  ```
  python3 determine.in-frame.stop-codon.py chm13.v1.0.pep.fasta chm13.draft_v1.0.gene_annotation.v4.gff3 | grep "chr22" > chm13.draft_v1.0.chr22.in-frame.stop-codon.bed
  ```

**2.** Identifying common variants from Chr22 VCF file (Aditi, Muhamad, Bryce, Tiancheng)

  -  Variant call filtering criteria (SNP, AAF > 5%)

  ```
  bcftools view -i "INFO/AF > 0.05" 1kgp.chr22.recalibrated.snp_indel.pass.withafinfo.vcf > 1kgp.chr22.recalibrated.snp_indel.pass.withafinfo.filtered_5%.vcf
  ```

**3.** Checking for overlaps between common variants from chr22 VCF file with stop codon sites identified from CHM13 fasta. Also, it will be checked if there are inconsistent nonsense variants or ORFs between CHM13 & hg38, which requirs RiboSeq validation (Anastasia, ChunHsuan)

  -  Picking up common variants

  ```
  (#$%$#%$%)
  ```

  -  Riboseq-validation (for the ORF & inframe-stop-codon sites where variants located)

  ```
  To download paired RNASeq.fastq and RiboSeq.fastq
  ```

  ```
  To quality control fastq files 
  (adapter trimming, remove low quality reads, etc.) 
  ```

  ```
  To align the reads to CHM13 refgenome 
  (get RNASeq.bam and RiboSeq.bam)
  ```

  ```
  To peak calling for read coverage of aligned bam files 
  (This will reveal real ORF and stop-codon sites by comparing RNAseq & Riboseq at same time.)
  ```

  ```
  To validate our targeted variant sites by the peak calling results, 
  and to clasify them into true ones and false ones.
  ```

#### IIII. Annotate common variants with ClinVar (Anastasia)



#### IV. Results:

**1.** Statistics of allele frequency

**2.** Postions of in-frame stop codons

    https://github.com/collaborativebioinformatics/popchrom/blob/main/data/chm13.draft_v1.0.chr22.in-frame.stop-codon.bed

**3.** Biologically annotated variants (CHM13 based).

**4.** Statistical visualization


## Required Data
- VCF files (vcf) 
  
  https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/by_chr/
  
- Reference sequence (CHM13 fasta & Hg38 fasta)
  
  https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/
  
- Gene annotation file (gff3, bed, others)

  https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/

- Example input: 

  Chr22.vcf - Variant called using CHM13 as reference
  
- Example output: 

  annotated variants.txt



## References

- Sergey Aganezov et.al, A complete reference genome improves analysis of human genetic variation. bioRxiv 2021.07.12.452063; doi: https://doi.org/10.1101/2021.07.12.452063
- T2T Consortium: https://sites.google.com/ucsc.edu/t2tworkinggroup/technology?authuser=0

## Acknoledgement

- Fritz J. Sedlazeck
- Ben Busby