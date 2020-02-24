# **Report 6: The Optimum K-mer and sim_threshold**

## *Summary*

1. [x] Abundance plot:
    1. References:
        - `extractedGTF_gencode.v33.transcripts.fa`
        - `longIso_extractedGTF_gencode.v33.transcripts.fa`
        - `nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa`
        - `noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa`
    2. k-mers to test: {25,31,41,51,75,81,91}
    3. Scatter plot with x-axis:kmer_size, y-axis:number of unique kmers with abundance > 1.
    4. By visual inspection, select the best k-mer size for the next steps.
2. [x] Create cDBG with the optimum kmer size from step #1.
3. [x] Clustering similarity threshold:
    1. [x] Run the CD-HIT-EST on the unitigs at different similarity threshold {91, 93, 95, 97, 99}.
    2. [x] Scatter plot with x-axis: similarity threshold, y-axis: number of clusters with size > 1.

## **1. Optimum k-mer size**

```bash

REF1=extractedGTF_gencode.v33.transcripts.fa
REF2=longIso_extractedGTF_gencode.v33.transcripts.fa
REF3=nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
REF4=noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa

for REF_FASTA in $REF1 $REF2 $REF3 $REF4
do
    HISTO_TSV=${REF_FASTA}_kmers_histogram.tsv
    touch ${HISTO_TSV}

    # 1.1 Jellyfish kmer count
    for K in 25 31 41 51 75 81 91;
    do
        jellyfish count -m ${K} -s 250M -t 10 -C ${REF_FASTA};
        UNIQ_KMERS=$(jellyfish dump mer_counts.jf  | awk -F">" '{if($2>1)a+=1}END{print a}');
        echo -e "${K}\t${UNIQ_KMERS}" >> ${HISTO_TSV};
        rm mer_counts.jf;
    done

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

```bash
python nested_kCount_histo.py
```

### 1.4.1 Plot 1 [Interactive](./plots/agg/plotly_histo.html)

![](./plots/agg/agg_kmers_histo.png?raw=true)

### 1.4.1 Plot 2 (Log scaled) [Interactive](./plots/agg/plotly_histo_log.html)

![](./plots/agg/agg_kmers_histo_log.png?raw=true)

---

## 2. Constructing cDBG for the reference files

```bash
OPTIMUM_KSIZE=75

REF1=extractedGTF_gencode.v33.transcripts.fa
REF2=longIso_extractedGTF_gencode.v33.transcripts.fa
REF3=nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
REF4=noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa

for REF_FASTA in $REF1 $REF2 $REF3 $REF4
do

    bcalm -kmer-size ${OPTIMUM_KSIZE} -max-memory 12000 -out cDBG_${REF_FASTA} -in ${REF_FASTA} &> log_cDBG_${REF_FASTA}
    rm *glue* *h5

done

```

## 3. CDHIT Clustering

```bash

REF1=extractedGTF_gencode.v33.transcripts.fa
REF2=longIso_extractedGTF_gencode.v33.transcripts.fa
REF3=nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
REF4=noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa

for REF_FASTA in $REF1 $REF2 $REF3 $REF4
do

    WORD_SIZE=9
    for SIM in 0.91 0.93 0.95
    do
        cd-hit-est -i cDBG_${REF_FASTA}.unitigs.fa -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_cDBG_${REF_FASTA} -d 0 -T 0 -M 12000 &> log_cdhit_${SIM}_${REF_FASTA}.log
    done

    WORD_SIZE=11
    for SIM in 0.97 0.99
    do
        cd-hit-est -i cDBG_${REF_FASTA}.unitigs.fa -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_cDBG_${REF_FASTA} -d 0 -T 0 -M 12000 &> log_cdhit_${SIM}_${REF_FASTA}.log
    done

done


```

## 3.2 Visualization of cd-hit clusters with size > 1 seq

```python
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

```

### 3.2.1 Plots and Data

```txt
notes:
    - Number of clusters = number of representative sequences from CDHIT
```

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
