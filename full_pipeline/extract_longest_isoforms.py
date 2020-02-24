import sys
from itertools import groupby
import textwrap
import os

"""
Extracting the longest isoforms from the CHR Basic annotation from gencode and the whole transcriptome file from Gencode.
"""

def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


fasta_file = sys.argv[1]
gene_to_transcripts = dict()

for header, seq in fasta_iter(fasta_file):
    header_info = header.split('|')
    transcript_id = header_info[0]
    gene_id = header_info[1]
    seq_len = len(seq)
    if gene_id not in gene_to_transcripts:
        gene_to_transcripts[gene_id] = dict()
    
    gene_to_transcripts[gene_id][transcript_id] = seq_len

longest_isoforms = set()
    
for gene, transcripts in gene_to_transcripts.items():
    longest_isoforms.add(max(transcripts, key=transcripts.get))

new_fasta = "longIso_" + os.path.basename(fasta_file)

with open(new_fasta, 'w') as fasta_writer:
    for header, seq in fasta_iter(fasta_file):
        transcript_id = header.split('|')[0]
        if transcript_id in longest_isoforms:
            fasta_writer.write(f">{header}\n")
            fasta_writer.write(textwrap.fill(seq, 60) + '\n')