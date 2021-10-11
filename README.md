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

After two decades of refinements, the human reference genome (GRCh38) has become more accurate and mostly complete. However, there are still hundreds of unresolved gaps persist, and no single chromosome has been finished from end to end because of the existence of highly repeated sequences (called transposable elements). Foutunatedly, by using the high-coverage & ultra-long-read technologies, several scientific groups have presented a new human genome assembly that surpasses the continuity of GRCh38, along with a gapless, telomere-to-telomere assembly of a human chromosome, which is called CHM13 reference genome. And for this brand new human reference genome, the precise annotation of variants and genes is required.

## Input:

- vcf files (using CHM13 as reference already).

## Outputs: 

- Correctly annotated variants sites.

## Methods 

I. Data Acquisition and Preprocessing:

1. Downloading the raw CHM13 fasta file.

2. Downloading the raw population based VCF (CHM13 fasta based).

3. Downloading the raw Annotation file (gff3).

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
<img width="1200" alt="flowchart" src="https://github.com/collaborativebioinformatics/popchrom/blob/main/others/flowchart_version2.png">
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

## Statistical Visualization of Example data

#### Input :

Example:<br/>
<img width="750" alt="examples......" src="https://github.com/collaborativebioinformatics/.png">

#### output :

Example:<br/>
<img width="750" alt="examples......" src="https://github.com/collaborativebioinformatics/.png">

## Appendix

## References 

- DNANexus documentation: https://documentation.dnanexus.com/developer/apps/execution-environment/connecting-to-jobs
- T2T Consortium: https://sites.google.com/ucsc.edu/t2tworkinggroup/technology?authuser=0
