# **Report 6: The Optimum K-mer**

## *Summary*

1. [x] Abundance plot:
    1. Reference: `longIso_extractedGTF_gencode.v33.transcripts.fa`
    2. k-mers to test: {25,31,41,51,75,81}
    3. Scatter plot with x-axis:kmer_size, y-axis:number of unique kmers with abundance > 1.
    4. By visual inspection, select the best k-mer size for the next steps.
2. [x] Create cDBG with the optimum kmer size from step #1.
3. [ ] Clustering similarity threshold:
    1. Run the CD-HIT-EST on the unitigs at different similarity threshold {90,93,96,99}.
    2. Scatter plot with x-axis: similarity threshold, y-axis: number of clusters with size > 1.

## **1. Optimum k-mer size**

```bash

# 1.1 Jellyfish count

HISTO_TSV=1_kmers_histogram.tsv
touch ${HISTO_TSV}
echo -e "kSize\tkmers>1\ttotal_kmers" >> ${HISTO_TSV}

for K in 25 31 41 51 75 81;
do
    jellyfish count -m ${K} -s 250M -t 10 -C longIso_extractedGTF_gencode.v33.transcripts.fa;
    TOTAL_KMERS=$(jellyfish dump mer_counts.jf  | grep ">" | wc -l);
    UNIQ_KMERS=$(jellyfish dump mer_counts.jf  | grep ">" | grep -o "[[:digit:]]*" | awk -F: '$NF+0 > 1' | wc -l);
    echo -e "${K}\t${UNIQ_KMERS}\t${TOTAL_KMERS}" >> ${HISTO_TSV};
    rm mer_counts.jf;
done


```

### 1.3 Histogram

```json
{'25': 6225004,
 '31': 5825143,
 '41': 5308103,
 '51': 4930689,
 '75': 4335645,
 '81': 4226369}
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

![](kmers_histo.png?raw=true)

### 1.4.1 Plot 2 (Log scaled)

![](kmers_histo_log.png?raw=true)

---

## 2. Constructing cDBG for the reference file longIso_extractedGTF_gencode.v33.transcripts.fa

```bash
OPTIMUM_KSIZE=75
bcalm -kmer-size ${OPTIMUM_KSIZE} -max-memory 12000 -out cDBG_longIso_extractedGTF_gencode.v33.transcripts -in longIso_extractedGTF_gencode.v33.transcripts.fa &> cDBGlongIso_extractedGTF_gencode.v33.transcripts.log

rm *glue* *h5

```

## 3. CDHIT Clustering

```bash

WORD_SIZE=9
for SIM in 0.90 0.93
do

    cd-hit-est -i cDBG_longIso_extractedGTF_gencode.v33.transcripts.unitigs.fa -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_cDBG_longIso_extractedGTF_gencode.v33.transcripts.unitigs -d 0 -T 0 -M 12000 &> cdhit_${SIM}.log

done

WORD_SIZE=11
for SIM in 0.96 0.99
do

    cd-hit-est -i cDBG_longIso_extractedGTF_gencode.v33.transcripts.unitigs.fa -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_cDBG_longIso_extractedGTF_gencode.v33.transcripts.unitigs -d 0 -T 0 -M 12000 &> cdhit_${SIM}.log

done


```