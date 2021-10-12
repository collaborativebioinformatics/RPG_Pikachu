#!/usr/bin/env python
import glob
import os
from snakemake.utils import min_version


shell.prefix("set -o pipefail; umask 002; ")  # set g+w


#Snakemake config
min_version("5.5")
configfile: "config.yaml"

if not os.path.exists("filtered_vcf"):
    os.makedirs("filtered_vcf")

##include rules (can be commented out)

include: "rules/filter_vcf.smk"
include: "rules/vcfconsensus.smk"
include: "rules/find_stopcodons.smk"


chr_list = list(range(1,23))+["X", "Y"]
chr_id= ["chr" + str(i) for i in chr_list]

rule all:
    input:
        "all_chr_fixed_header.vcf",
        expand("input_vcf/{chromosome}.vcf", chromosome = chr_id),
        expand("filtered_vcf/{chromosome}_filtered.vcf", chromosome = chr_id),
        "filtered_vcf/common_annotated.vcf.gz",
        "output/CHM13_new.fasta",
        "output/chm13.new.in-frame.stop-codon.bed"

rule fix_header:
    input:
        vcf = config["vcf"]
    output:
        "all_chr_fixed_header.vcf"
    conda:
        "envs/utils.ymal"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        gatk FixVcfHeader \
          I={input} \
          O={output}
        """


rule split_by_chr:
    input:
        vcf = "all_chr_fixed_header.vcf"
    output:
        "input_vcf/{chromosome}.vcf"
    params:
        chr = "{chromosome}"
    conda:
        "envs/utils.ymal"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        gatk SelectVariants \
          -R {config[ref_genome]} \
          -V {input.vcf} \
          -L {params.chr} \
          --exclude-filtered true \
          -O {output}
        """

rule mergevcf:
    input:
        expand("filtered_vcf/{chromosome}_filtered.vcf.gz", chromosome=chr_id)
    output:
        "filtered_vcf/common.vcf.gz"
    conda:
        "envs/utils.ymal"
    threads: 10
    resources:
        nodes = 1
    shell:
        """
        bcftools concat {input} \
                 --threads {threads} \
                 -O z \
                 -o {output} && tabix -p vcf {output}       
        """

rule annotatevcf:
    input:
        vcf = "filtered_vcf/common.vcf.gz",
        vcf_clinvar = {config["vcf_clinvar"]}
    output:
        "filtered_vcf/common_annotated.vcf.gz"
    conda:
        "envs/utils.ymal"
    threads: 10
    resources:
        nodes = 1
    shell:
        """
        gatk CreateSequenceDictionary -R {config[ref_genome]}
        gatk VariantAnnotator \
          -R {config[ref_genome]} \
          -V {input.vcf} \
          -O {output} \
          -resource:clinvar {input.vcf_clinvar} -E clinvar.CLNSIG
        """

