# *Report7: Extension by SSAKE*

## **Summary**

1. [x] Create cDBG with the optimum kmerSize from the raw reads.
2. [x] CD-HIT-EST the #1 cDBG with the **optimum** similarity threshold.
3. [*] Extend the representative sequences from #2 using [SSAKE](https://github.com/bcgsc/SSAKE).
4. [*] Get the longest from SSAKE write a new fasta file called ssake_after_75.
5. [-] Compare the #3 fasta with the reprs from CDHIT clusters: Length of the two sequences and report the stats.
6. [-] Create cDBG ssake_after_75 for the SSAKE reps.
7. [ ] Rerun #2 to 4 for one extra time, rename the output fasta into ssake_after_75_2.

## **2 CDHIT Clustering  SRR11015356_k75.unitigs with sim=0.99**

### **2.1 Clustering**

```bash
BEFORE_75_UNITIGS=/media/mabuelanin/26c644c2-7c3c-4a3d-8172-596cf2932040/home/mabuelanin/symbolic/data/drtamer_data/SRR11015356_k75.unitigs.fa
cd-hit-est -i ${BEFORE_75_UNITIGS} -n 11 -c 0.99 -o cdhit_0.99_SRR11015356_k75.unitigs -d 0 -T 0 -M 12000
```

### **2.2 Results**

11,824,622 unitigs clustered into 11,019,543 clusters.
> Total CPU time 21253.84
---

## **3 Extending the clusters' representative sequences by SSAKE**

### **3.1 Preprocessing to output lists of IDs**

- The following python script will generate multiple files:
    1. **ssake_target_seqs.txt:** contains pipe separated for ssake target sequences.
    2. **ssake_all_seeds.txt:** contains list of seed sequences in the same order of file 1.
    3. **all_exluded_seqs.txt:** contains union list of seeds and targets to be extracted from the unitigs fasta file.
    4. **passed_reprs.txt:** contains a list of representative sequence IDs of single-sequence-clusters.

- `paste ssake_all_seeds.txt ssake_target_seqs.txt | head`:

```tsv
>10171743	>9410083|>10832746
>8633421	>8119486|>10377219
>10087844	>8120313|>11365534
>7309652	>8553822|>10816700
>8454780	>2085730|>10155878
>7920550	>8072610|>9032344
```

---

```python
from itertools import groupby
import re
import os


def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

cdhitClustersFile = "cdhit_0.99_SRR11015356_k75.unitigs.clstr"


large_clusters = dict()
single_IDs = []
large_IDs = set()

print("Parsing ...")
for clusterID, sequences in fasta_iter(cdhitClustersFile):
    if sequences.count('>') > 2:
        clusterID = clusterID.split()[-1]
        reprs = re.findall(r'>(\d+)\.{3}\s\*', sequences)[0]
        others = re.findall(r'>(\d+)\.{3}\sat', sequences)
        large_clusters[clusterID] = [reprs, others]
        large_IDs.add(reprs)
        for _ in others:
            large_IDs.add(_)

    else:
        reprs = re.findall(r'>(\d+)\.{3}\s\*', sequences)[0]
        single_IDs.append(f"{reprs}")
        #single_IDs.append(clusterID.split()[-1])

print("Writing ...")

with open("passed_reprs.txt", 'w') as passedRepsWriter:
    for _id in single_IDs:
        passedRepsWriter.write(f">{_id}\n")

with open("ssake_all_seeds.txt" , 'w') as ssakeAllSeedsWriter, open('ssake_target_seqs.txt', 'w') as ssakeTargetsWriter, open("all_exluded_seqs.txt", 'w') as allExcludedWriter:
    for clusterID, sequences in large_clusters.items():

        ssakeAllSeedsWriter.write(f">{sequences[0]}\n")

        _line = [f">{sq}" for sq in sequences[1]]
        ssakeTargetsWriter.write('|'.join(_line) + '\n')

        allExcludedWriter.write(f">{sequences[0]}\n")
        for s in sequences[1]:
            allExcludedWriter.write(f">{s}\n")
```

### **3.2 preparing ssake files**

```bash

# Extracting all unready sequences from large clusters from the whole unitigs file
cat SRR11015356_k75.unitigs.fa | grep -A1 --no-group-separator -E -w -f all_exluded_seqs.txt > largeClusters_SRR11015356_k75.unitigs.fa

# Extract all the passed representatives of single-sequence clusters
cat SRR11015356_k75.unitigs.fa | grep -A1 --no-group-separator -w -f passed_reprs.txt > passedReprs_SRR11015356_k75.unitigs.fa

# Extract original seeds from the whole unitigs file to compare the difference later on.
cat SRR11015356_k75.unitigs.fa | grep -A1 --no-group-separator -w -f ssake_all_seeds.txt > originalSeeds_SRR11015356_k75.unitigs.fa


```

## **2.2 SSAKE Run**

```bash

LC_ALL=C

OUTPUT_REPRS="ssake_repr_SRR11015356_k75.unitigs.fa"

rm -rf ${OUTPUT_REPRS}
touch ${OUTPUT_REPRS}

paste ssake_target_seqs.txt ssake_all_seeds.txt | while read targets seed;
do
  cat largeClusters_SRR11015356_k75.unitigs.fa | fgrep -A1 -w ${seed} > tmp.seed.fa;
  cat largeClusters_SRR11015356_k75.unitigs.fa | grep -A1 --no-group-separator -E -w ${targets} > tmp.targets.fa;
  ./ssake/SSAKE -f tmp.targets.fa -w 1 -s tmp.seed.fa -i 0 -o 1 -r 0.6 -y 1 -h 1 -q 1 -b ssake_contigs -v 0 &> /dev/null;
  cat ssake_contigs_singlets.fa >> ${OUTPUT_REPRS};
done

rm ssake_contigs*
rm -rf tmp.seed.fa tmp.targets.fa

```

---

## **4. Summary stats of before & after SSAKE**

```bash
seqkit stats ssake_repr_SRR11015356_k75.unitigs.fa originalSeeds_SRR11015356_k75.unitigs.fa
```

```tsv
waiting for the output ...
```

---

## **5. Create cDBG ssake_after_75 for the SSAKE reps**

### 5.1 merging ssake_seads with the passed cdhit seeds **curated_cdhit_0.99_SRR11015356_k75.unitigs.fa**

```bash
cat ssake_repr_SRR11015356_k75.unitigs.fa passedReprs_SRR11015356_k75.unitigs.fa > curated_cdhit_0.99_SRR11015356_k75.unitigs.fa
```

### 5.2 Creating the cDBG

```bash
bcalm -kmer-size 75 -max-memory 12000 -out ssake_before75_after75 -in curated_cdhit_0.99_SRR11015356_k75.unitigs.fa &> curated_cdhit_0.99_SRR11015356_k75_k75.log
```
