#!/usr/bin/env python
import glob
from pathlib import Path


rule filter_vcf:
    input:
        vcf = "final_vcfs/{chromosome}.vcf", 
    params:
        chr = "{chromosome}"
    output:
        "filtered_vcf/{chromosome}_filtered.vcf"
    shell:
        """
        bcftools view -i 'INFO/AF > 0.05' {input.vcf} > {input.vcf}_filteredAF.vcf
        bcftools view --max-alleles 2 --exclude-types indels {input.vcf}_filteredAF.vcf > {output}
        rm {input.vcf}_filteredAF.vcf
        """

rule filter_vcf_index:
    input:
        vcf = "final_vcfs/{chromosome}_filtered.vcf", 
    params:
        chr = "{chromosome}"
    output:
        "filtered_vcf/{chromosome}_filtered.vcf.gz"
    shell:
        """
        bgzip -c {input.vcf} > {output}
        tabix -p vcf {output}
        """







