#!/usr/bin/python3

import re, sys
from Bio import SeqIO

def read_gff(infile):
    f = open(infile, 'r')
    gff_dict = dict()
    trans_dict = dict()
    
    for line in f:
        line = line.strip()
        if line.startswith('#'): continue
        content = line.split()
        if content[2] == "transcript":
            gid = re.search(r'source_gene=(.*?);', content[8]).group(1)
            tid = re.search(r'transcript_id=(.*?);', content[8]).group(1)
            gff_dict[tid] = [content[0], content[3]]
            trans_dict[tid] = gid
            
    return gff_dict, trans_dict

def run(fastafile, gff_file):
    gff_dict, trans_dict = read_gff(gff_file)
     
    f = open(fastafile, 'r')
    for record in SeqIO.parse(f, 'fasta'):
        seqid = record.id
        sequence = record.seq

        for index, base in enumerate(sequence):
            if base == ".":
                relative_start = index * 3
                start_pos = int(gff_dict[seqid][1]) + relative_start - 1
                
                print(gff_dict[seqid][0] + "\t" + str(start_pos) + "\t" + str(start_pos + 3) + "\t" + seqid + ":" + trans_dict[seqid])
    f.close()

if __name__ == "__main__":
    protein_file = sys.argv[1]
    gff_file = sys.argv[2]
    run(protein_file, gff_file)
