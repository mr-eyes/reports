import sys
from itertools import groupby
import textwrap
import os

"""
This works with the CHR Basic annotation from gencode and the whole transcriptome file from Gencode.
"""

def fasta_iter(fasta_name):
    """
    Thanks to https://www.biostars.org/p/710
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name, 'r')
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


if len(sys.argv) < 3:
    sys.exit("run python extract_transcripts_from_gtf.py <gtf> <fasta>")

gtf_file = sys.argv[1]
fasta_file = sys.argv[2]

transcripts = list()

with open(gtf_file, 'r') as gtf_reader:
    for line in gtf_reader:
        if line[0] == '#': continue
        part1, part2 = tuple(line.split("gene_id"))
        record_type = part1.split('\t')[2]
        if record_type == 'transcript':
            transcript_id = part2.split(";")[1].split(' ')[2].replace('"','')
            transcripts.append(transcript_id)

new_fasta = "extractedGTF_" + os.path.basename(fasta_file)

with open(new_fasta, 'w') as fasta_writer:
    for header, seq in fasta_iter(fasta_file):
        transcript_id = header.split('|')[0]
        if transcript_id in transcripts:
            fasta_writer.write(f">{header}\n")
            fasta_writer.write(textwrap.fill(seq, 60) + '\n')