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

### 1.3 Histogram visualization

<h3><b>The following plots represents the number of kmers with abundance > 1 </b></h3>

```bash
python plotly_nested_kCount_histo.py
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

```bash
python plotly_clusters_size.py
```

### 3.3 Plots


### 3.3.1 Plot 1 [Interactive](./plots/agg/plotly_cdhit.html)

![](./plots/agg/plotly_cdhit.png?raw=true)


### 3.3.2 Plot 2 (Log scaled) [Interactive](./plots/agg/plotly_cdhit_log.html)

![](./plots/agg/plotly_cdhit_log.png?raw=true)

### 3.3.3 Plot 3 (total clusters) [Interactive](./plots/agg/plotly_cdhit_total.html)

![](./plots/agg/plotly_cdhit_total.png?raw=true)
