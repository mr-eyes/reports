"""
Annotating swissProt sequences by the m8 MMseqs2 output file                                      #

[0]query [1]target [2]pident [3]alnlen [4]mismatch [5]gapopen [6]qstart [7]qend [8]tstart [9]tend [10]evalue [11]bits 
P62807	Q6ZWY9	1.000	126	0	0	1	126	1	126	3.978E-71	236
P62807	Q5R893	1.000	126	0	0	1	126	1	126	3.978E-71	236
"""

from itertools import groupby
import textwrap
import sys

if len(sys.argv) < 2:
    sys.exit("run: python transform_m8.py <m8_file> <fasta_file>")

m8_file = sys.argv[1]
#fasta_file = sys.argv[2]

species_to_human = dict()
species_to_AG = dict()
pair_to_pident = dict()

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

# Parsing
with open(m8_file, 'r') as m8:
    for line in m8:
        line = line.strip().split()
        human = line[1]
        species = line[0]
        evalue = float(line[10])
        pident = float(line[2])
        pair_to_pident[(human,species)] = pident
        if species not in species_to_human:
            species_to_human[species] = dict()
                    
        species_to_human[species][human] = evalue


for species, human_scores in species_to_human.items():
    species_to_AG[species] = min(human_scores, key=human_scores.get)

xx = list()

for species, AG in species_to_AG.items():
    xx.append((species, AG))
#    print(f"{species}\t{AG}\n")


with open(m8_file, 'r') as m8:
    for line in m8:
        l = line
        line = line.strip().split()
        human = line[1]
        species = line[0]
        if (species, human) in xx:
            print(l, end='')

exit()
# Rewriting
with open(f"annotated_{fasta_file}", 'w') as new_fasta:
    for header, seq in fasta_iter(fasta_file):
        species_seqID = header.split('|')[1]
        best_humanGene = species_to_AG.get(species_seqID, "NaN")
        header += f"AG={best_humanGene}"
        new_fasta(textwrap.fill(seq, 60))
        