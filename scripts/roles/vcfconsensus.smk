#!/usr/bin/env python
import glob
from pathlib import Path


rule vcf_consensus:
    input:
        ref = config["ref_genome"],
        vcf = "filtered_vcf/common.vcf.gz", 
    output:
        "output/CHM13_new.fasta"
    shell:
        """
        cat {input.ref} | vcf-consensus {input.vcf} > {output}
        """


