# **Report 9: kProcessor test**

- [x] Create cDBG @ k=75 for raw reads
- [x] Index the cDBG  @ k=25 using MQF **[FAILED BY EVERY WAY]**
- [ ] Try indexing a small sample


## 1. Indexing the full unitigs file with k=25

### 1.1 Generate names file

```bash
grep ">" SRR11015356_k75.unitigs.fa | cut -c2- |  awk -F' ' '{print $0"\t"$1}' > SRR11015356_k75.unitigs.fa.names
``  

### 1.2 Indexing

# Failed and reported in [<h3>kProcessor:#38</h3>](https://github.com/dib-lab/kProcessor/issues/38)

--

## 2. Indexing representative sequences only of the unitigs

### 2.1 Generate names file

```bash
grep ">" reps_unitigs_SRR11015356_k75.fa | cut -c2- |  awk -F' ' '{print $0"\t"$1}' > reps_unitigs_SRR11015356_k75.fa.names
```

### 2.2 Indexing

```python
import kProcessor as kp
fasta_file="/home/mabuelanin/Desktop/kexpression_experiment/symbolic/data/drtamer_data/reps_unitigs_SRR11015356_k75.fa"
names_file=fasta_file + ".names"


kSize = 25
Q = 30
mode = 1 # kmers
chunk_size = int(1e4)

kf = kp.kDataFrameMQF(kSize, Q, mode)

print("Indexing ...")

ckf = kp.index(kf, {"mode" : mode}, fasta_file , chunk_size, names_file)

print("Serializing the index ...")

ckf.save("idx_k25_reps_unitigs_SRR11015356_k75")
```