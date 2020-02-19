from itertools import groupby
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd


def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq



histo = {
    "large_clusters" : dict(),
    "total_clusters" : dict(),
}


for file in sorted(glob("./*.clstr")):
    sim_threshold = file.split("_")[1]
    total_clusters = 0
    large_clusters = 0

    for header, seq in fasta_iter(file):
        size = seq.count('nt')
        total_clusters += 1
        if size > 1:
            large_clusters += 1

    histo["total_clusters"][sim_threshold] = total_clusters
    histo["large_clusters"][sim_threshold] = large_clusters


pd.DataFrame(histo).plot(kind='bar')
plt.show()