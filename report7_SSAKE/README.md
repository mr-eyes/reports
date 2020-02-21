# *Report7: Extension by SSAKE*

## **Summary**

1. [x] Create cDBG with the optimum kmerSize from the raw reads.
2. [x] CD-HIT-EST the #1 cDBG with the **optimum** similarity threshold.
2. [ ] Extend the representative sequences from #2 using [SSAKE](https://github.com/bcgsc/SSAKE).
3. [ ] Get the longest from SSAKE write a new fasta file called ssake_after_75.
4. [ ] Compare the #3 fasta with the reprs from CDHIT clusters: Length of the two sequences and report the stats.
4. [ ] Create cDBG ssake_after_75 for the SSAKE reps.
5. [ ] Rerun #2 to 4 for one extra time, rename the output fasta into ssake_after_75_2.

## **2 CDHIT Clustering  SRR11015356_k75.unitigs with sim=0.99**

### **2.1 Clustering**

```bash
BEFORE_75_UNITIGS=/media/mabuelanin/26c644c2-7c3c-4a3d-8172-596cf2932040/home/mabuelanin/symbolic/data/drtamer_data/SRR11015356_k75.unitigs.fa
cd-hit-est -i ${BEFORE_75_UNITIGS} -n 11 -c 0.99 -o cdhit_0.99_SRR11015356_k75.unitigs -d 0 -T 0 -M 12000
```

### **2.2 Results**

11,824,622 unitigs clustered into 11,019,543 clusters.

---

## **3 Extending the clusters' representative sequences by SSAKE**

```bash

SSAKE -f all_except_representitive.fa -w 1 -s representitive.fa -i 0 -o 1 -r 0.6 -y 1 -h 1 -q 1


SSAKE -f all_except_representitive.fa \
-w 1 \
-s representitive.fa -i 0 \
-o 1 \
-r 0.6 \  ## I am not sure (default is 0.7: Minimum base ratio used to accept a overhang consensus base. Higher -r value lead to more accurate contig extension.
-y 1 \
-h 1 \
-q 1

```