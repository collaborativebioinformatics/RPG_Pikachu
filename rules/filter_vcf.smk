#!/usr/bin/env python
import glob
from pathlib import Path

shell.prefix("set -o pipefail; umask 002; ")  # set g+w
rule filter_vcf:
    input:
        vcf = "input_vcf/{chromosome}.vcf"
    output:
        "filtered_vcf/{chromosome}_filtered.vcf"
    conda:
        "../envs/utils.ymal"
    shell:
        """
        bcftools view -i 'INFO/{config[af_flag]} > 0.05' {input.vcf} > {input.vcf}_filteredAF.vcf
        bcftools view --max-alleles 2 --exclude-types indels {input.vcf}_filteredAF.vcf > {output}
        rm {input.vcf}_filteredAF.vcf
        """

rule filter_vcf_index:
    input:
        vcf = "filtered_vcf/{chromosome}_filtered.vcf", 
    output:
        "filtered_vcf/{chromosome}_filtered.vcf.gz"
    conda:
        "../envs/utils.ymal"
    shell:
        """
        bgzip -c {input.vcf} > {output}
        tabix -p vcf {output}
        """







