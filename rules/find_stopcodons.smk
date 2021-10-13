#!/usr/bin/env python
import glob
from pathlib import Path

shell.prefix("set -o pipefail; umask 002; ")  # set g+w
rule get_inframe_stopCodon:
    input:
        ref=config["ref_genome"], 
        gff=config["gff"]
    conda:
        "../envs/utils.ymal"
    output:
        "output/chm13.stop-codon.bed"
    shell:
        """
        gffread -g {input.ref} {input.gff} -y pep.fasta
        python3 determine.in-frame.stop-codon.py pep.fasta {input.gff} > {output}
        rm new.pep.fasta
        """

