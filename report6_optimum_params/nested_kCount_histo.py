import matplotlib.pyplot as plt
import pandas as pd
import os


kmers = [31,41,51,75,81]
data = dict()
for k in kmers:
    data[f"k{k}"] = dict()


ref1="test_data/ref1.tsv"
ref2="test_data/ref2.tsv"

for ref in [ref1, ref2]:
    with open(ref , 'r') as tsv_reader:
        for line in tsv_reader:
            line = line.strip().split()
            kSize = line[0]
            uniq = int(line[1])
            data[f"k{kSize}"][os.path.basename(ref)] = uniq



pd.DataFrame(data).plot(kind='bar')
plt.show()