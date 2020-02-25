# **Report 9: kProcessor test**



- [x] Create cDBG @ k=75 for raw reads
- [x] Index the cDBG  @ k=25 using MQF **[FAILED BY EVERY WAY]**


## Generate names file

```bash
cat SRR11015356_k75.unitigs.fa | grep ">" | awk -F\, '{print (substr($0,2))"\t"substr($0,2,1)}' > SRR11015356_k75.unitigs.fa.names
``  

```python
import kProcessor as kp

fasta_file="SRR11015356_k75.unitigs.fa"
names_file= fasta_file + ".names"


kSize = 25

# Tried with 29,30,31,32,33,34
Q = 33

# kmers mode
mode = 1

chunk_size = int(1e4)

# Tried and failed
# kf = kp.kDataFrameMQF(kSize, Q, mode)

# Tried & failed immediately with small dataset, running forever with the large dataset.
# kf = kp.kDataFrameBMQF(kSize)

# Failed but the problem (stack trace) was in the data structure holding the colors.
kf = kp.kDataFramePHMAP(kSize)

print("Indexing ...")

ckf = kp.index(kf, {"mode" : mode}, fasta_file , chunk_size, names_file)

print("Serializing the index ...")

ckf.save("idx_k25_SRR11015356_k75.unitigs")

```
