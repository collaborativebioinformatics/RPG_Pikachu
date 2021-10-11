# Telomere-to-telomere common alleles & abnormally avoided stop codons

## Contributors 
-  Elbay Aliyev - Lead 
-  Muhammad Sohail Raza - Writer
-  Shangzhe Zhang, Anastasia Illarionova, Prashant Ranjan, Tiancheng Xu, Aditi Sammi -　Tech support 
-  Bryce Kille, ChunHsuan LO - SysAdmin

## Goal 
- Identify rare variants that has stop-codons in CHM13 and replacing with common variant? Checking the fasta file directly.
- Identify stop codons that disagree with Ribo-seq analysis.

## Introduction 

## Input:

- vcf files lifted up from hg38 to CHM13.

## Outputs: 

- Correctly annotated variants site.

## Methods 

I. Data Acquisition and Preprocessing:

1. Download the row VCF (hg38 fasta based)

II. Core tasks:

1. Identifying stop codon sites in CHM13 fasta file (Shangzhe, Muhamad, Bryce)

2. Identifying common variants from Chr22 VCF file (Aditi, Muhamad, Bryce)

3. Checking for overlaps between common variants from chr22 VCF file with stop codon sites identified from CHM13 fasta. Also, it will be checked if there are inconsistent nonsense variants between CHM13 & hg38, which requirs RiboSeq validation (Anastasia, ChunHsuan)

4. Annotate common variants with ClinVar (Anastasia, Shangzhe)

5. Flowchart creation (Anastasia, Muhamad)

6. Biological significance of the replacement selection/disease-related candidates. 

III. Outcome:

1. Biologically annotated variants.

## Installation 
Please use the DNAnexus workflows to use this tool. 

## Flowchart
<img width="1200" alt="flowchart" src="https://github.com/collaborativebioinformatics/.jpg">
(for the working pipelines)

## Required Data
- VCF files (vcf) 
  
- Reference sequence (fasta)

  CHM13 Fasta file, Hg38 Fasta file
  
- Annotated reference file (gff3)

  https://s3-us-west-2.amazonaws.com/human-pangenomi

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
