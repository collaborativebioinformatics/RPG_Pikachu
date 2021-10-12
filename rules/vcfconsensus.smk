#!/usr/bin/env python
import glob
from pathlib import Path

shell.prefix("set -o pipefail; umask 002; ")  # set g+w
rule vcf_consensus:
    input:
        ref = config["ref_genome"],
        vcf = "filtered_vcf/common.vcf.gz"
    output:
        "output/CHM13_new.fasta"
    conda:
        "../envs/utils.ymal"
    shell:
        """
        cat {input.ref} | vcf-consensus {input.vcf} > {output}
        """
