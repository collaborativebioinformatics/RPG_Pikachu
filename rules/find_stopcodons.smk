#!/usr/bin/env python
import glob
from pathlib import Path

shell.prefix("set -o pipefail; umask 002; ")  # set g+w
rule get_inframe_stopCodon:
    input:
        ref='output/CHM13_new.fasta', 
        gff=config["gff"]
    conda:
        "../envs/utils.ymal"
    output:
        "output/chm13.new.in-frame.stop-codon.bed"
    shell:
        """
        gffread -g {input.ref} {input.gff} -y > new.pep.fasta
        python3 determine.in-frame.stop-codon.py new.pep.fasta {input.gff} > {output}
        rm new.pep.fasta
        """

