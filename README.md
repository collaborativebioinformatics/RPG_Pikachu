# Telomere-to-telomere common alleles & abnormally avoided stop codons

## Contributors 
-  Elbay Aliyev - Lead 
-  Muhammad Sohail Raza - Writer
-  Shangzhe Zhang, Anastasia Illarionova, Prashant Ranjan, Tiancheng Xu, Aditi Sammi -　Tech support 
-  Bryce Kille, ChunHsuan LO - SysAdmin

## Goal 
- To identify rare variants that has stop-codons in CHM13 and recorrecting some of them as common variants. By Checking the fasta file directly.
- To identify stop codons that disagree with Ribo-seq analysis for validting the annotation.

## Introduction 

## Input:

- vcf files (using CHM13 as reference already).

## Outputs: 

- Correctly annotated variants sites.

## Methods 

I. Data Acquisition and Preprocessing:

1. Downloading the raw CHM13 fasta file.

2. Downloading the raw population based VCF (CHM13 fasta based).

3. Downloading the raw Annotation file (gff3).
4. 
II. Core tasks:

1. Identifying stop codon sites in CHM13 fasta file (Shangzhe, Muhamad, Bryce)

2. Identifying common variants from Chr22 VCF file (Aditi, Muhamad, Bryce)

3. Checking for overlaps between common variants from chr22 VCF file with stop codon sites identified from CHM13 fasta. Also, it will be checked if there are inconsistent nonsense variants between CHM13 & hg38, which requirs RiboSeq validation (Anastasia, ChunHsuan)

4. Annotate common variants with ClinVar (Anastasia, Shangzhe)

5. Flowchart creation (Anastasia, Muhamad)

6. Biological significance of the replacement selection/disease-related candidates. 

7. Statistical visualization.

III. Outcome:

1. Biologically annotated variants (CHM13 based).

## Installation 
Please use the DNAnexus workflows to use this tool. 

## Flowchart
<img width="1200" alt="flowchart" src="https://github.com/collaborativebioinformatics/popchrom/blob/main/others/flowchart_version1.png">
(for the working pipelines)

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

### Appendix

#### Input :

Example:<br/>
<img width="750" alt="PRS_value of each feature for the example person." src="https://github.com/collaborativebioinformatics/.png">

## References 

- DNANexus documentation: https://documentation.dnanexus.com/developer/apps/execution-environment/connecting-to-jobs
