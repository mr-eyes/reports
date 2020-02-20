# **Report 6: The Optimum K-mer**

## *Summary*

1. [x] Abundance plot:
    1. Reference: `nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa`
    2. k-mers to test: {25,31,41,51,75,81}
    3. Scatter plot with x-axis:kmer_size, y-axis:number of unique kmers with abundance > 1.
    4. By visual inspection, select the best k-mer size for the next steps.
2. [x] Create cDBG with the optimum kmer size from step #1.
3. [x] Clustering similarity threshold:
    1. [x] Run the CD-HIT-EST on the unitigs at different similarity threshold {90,93,96,99}.
    2. [x] Scatter plot with x-axis: similarity threshold, y-axis: number of clusters with size > 1.

## **1. Optimum k-mer size**

```bash

# 1.1 Jellyfish count

HISTO_TSV=1_kmers_histogram.tsv
touch ${HISTO_TSV}
echo -e "kSize\tkmers>1\ttotal_kmers" >> ${HISTO_TSV}

for K in 25 31 41 51 75 81;
do
    jellyfish count -m ${K} -s 250M -t 10 -C nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa;
    TOTAL_KMERS=$(jellyfish dump mer_counts.jf  | grep ">" | wc -l);
    UNIQ_KMERS=$(jellyfish dump mer_counts.jf  | awk -F">" '{if($2>1)a+=1}END{print a}');
    echo -e "${K}\t${UNIQ_KMERS}\t${TOTAL_KMERS}" >> ${HISTO_TSV};
    rm mer_counts.jf;
done


```

### 1.3 Histogram

```tsv
kSize	kmers>1	total_kmers
25	3466029	89464767
31	3139993	90421021
41	2739366	91196188
51	2463775	91422034
75	2066541	91080253
81	1999465	90920059
```

### 1.4 Histogram visualization

```python
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

```

### 1.4.1 Plot 1

![](run2_kmers_histo.png?raw=true)

### 1.4.1 Plot 2 (Log scaled)

![](run2_kmers_histo_log.png?raw=true)

---

## 2. Constructing cDBG for the reference file nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa

```bash
OPTIMUM_KSIZE=75
bcalm -kmer-size ${OPTIMUM_KSIZE} -max-memory 12000 -out cDBG_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts -in nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa &> cDBGlongIso_nonoverlap_extractedGTF_gencode.v33.transcripts.log

rm *glue* *h5

```

## 3. CDHIT Clustering

```bash

WORD_SIZE=9
for SIM in 0.90 0.93
do

    cd-hit-est -i cDBG_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.unitigs.fa -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_cDBG_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.unitigs -d 0 -T 0 -M 12000 &> cdhit_${SIM}.log

done

WORD_SIZE=11
for SIM in 0.96 0.99
do

    cd-hit-est -i cDBG_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.unitigs.fa -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_cDBG_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.unitigs -d 0 -T 0 -M 12000 &> cdhit_${SIM}.log

done


```

## 3.2 Visualization of cd-hit clusters with size > 1 seq

```python

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
    "total_clusters x10" : dict(),
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

    histo["total_clusters x10"][sim_threshold] = total_clusters // 10
    histo["large_clusters"][sim_threshold] = large_clusters


pd.DataFrame(histo).plot(kind='bar')
plt.show()

```

### 3.2.1 Plots and Data

```json
{
   "large_clusters":{
      "0.90":1651,
      "0.93":1626,
      "0.96":1624,
      "0.99":494
   },
   "total_clusters":{
      "0.90":11590,
      "0.93":12035,
      "0.96":12821,
      "0.99":15773
   }
}
```

![](run2_cdhit_histo_raw.png?raw=true)

![](run2_cdhit_histo_x10.png?raw=true)
