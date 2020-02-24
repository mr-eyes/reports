from itertools import groupby
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import os

def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq



REF1="extractedGTF_gencode.v33.transcripts.fa"
REF2="longIso_extractedGTF_gencode.v33.transcripts.fa"
REF3="nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa"
REF4="noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa"
REFS = [REF1, REF2, REF3, REF4]


histo = dict()
sims = [91,93,95,97,99]

for sim in sims:
    histo[f"{sim}%"] = dict()


for REF_FASTA in REFS:
    for sim in sims:
        file = f"clusters_0.{sim}_cDBG_{REF_FASTA}.clstr"
        large_clusters = 0

        for header, seq in fasta_iter(file):
            size = seq.count('nt')
            total_clusters += 1
            if size > 1:
                large_clusters += 1

        title = REF_FASTA.split('.')[0].split('_')[:-1]
        histo[f"{sim}%"][title] = large_clusters


pd.DataFrame(histo).plot(kind='bar')
plt.show()