import matplotlib.pyplot as plt
import pandas as pd
import os


kmers = [25, 31,41,51,75,81, 91]
data = dict()
for k in kmers:
    data[f"k{k}"] = dict()


refs = [
"extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv",
"longIso_extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv",
"nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv",
"noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv"
]

for ref in refs:
    with open(ref , 'r') as tsv_reader:
        for line in tsv_reader:
            line = line.strip().split()
            kSize = line[0]
            uniq = int(line[1])
            title = os.path.basename(ref).split('_')[0] # .split('_')[:-1]
            data[f"k{kSize}"][title] = uniq


pd.DataFrame(data).plot(kind='bar')
plt.show()