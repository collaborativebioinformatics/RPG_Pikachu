#!/usr/bin/env python
import glob
from pathlib import Path


rule get_pep_sequence:
    input:
        ref='output/CHM13_new.fasta', 
        gff=config["gff"]
    conda:
        "envs/utils.ymal"
    output:
        "output/chm13.new.pep.fasta"
    shell:
        "gffread -g {input.ref} {input.gff} -y > {output}

rule get_inframe_stopCodon:
    input:
        fa="output/chm13.new.pep.fasta",
        gff="data/chm13.draft_v1.0.gene_annotation.v4.gff3"
    conda:
        "envs/utils.ymal"
    output:
        "chm13.new.in-frame.stop-codon.bed"
    shell:
        "python3 determine.in-frame.stop-codon.py {input.fa} {input.gff} > {output}"
