import matplotlib.pyplot as plt
import pandas as pd

data = {"uniq_kmers_count" : dict()}

with open("1_kmers_histogram.tsv" , 'r') as tsv_reader:
    next(tsv_reader)
    for line in tsv_reader:
        line = line.strip().split()
        kSize = line[0]
        uniq = int(line[1])
        total = int(line[2])
        data["uniq_kmers_count"][kSize] = uniq


pd.DataFrame(data).plot(kind='bar')
plt.show()