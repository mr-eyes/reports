# Full pipeline

## Summary

## **1. Data preparation**

### 1.1 Downloading & Dumping

```bash

#1 Human Transcriptome
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz
gunzip gencode.v33.transcripts.fa.gz

#2 Basic gene annotation CHR GTF
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.basic.annotation.gtf.gz
gunzip gencode.v33.basic.annotation.gtf.gz

#3 Representitive sequences extractions
python extract_transcripts_from_gtf.py gencode.v33.basic.annotation.gtf gencode.v33.transcripts.fa

#4 Sample reads downloading & dumping
wget -c https://sra-download.ncbi.nlm.nih.gov/traces/sra51/SRR/010757/SRR11015356 -O SRR11015356.sra
fastq-dump --fasta 0 --split-files SRR11015356.sra

```

### 1.2 cDBG of "SRR11015356_1.fasta" and "SRR11015356_2.fasta" @ k=75

```bash
ls -1 *fasta > list_reads
bcalm -kmer-size 75 -max-memory 12000 -out SRR11015356_before_k75 -in list_reads
```

### 1.3 Clustering the cDBG at similarity threshold = 0.95

```bash
cd-hit-est -i SRR11015356_before_k75.unitigs.fa -n 11 -c 0.95 -o clusters_SRR11015356_before_k75 -d 0 -T 0 -M 12000
```

### 1.4 Exporting representative sequences only from the unitigs.fa files

```bash
cat clusters_SRR11015356_before_k75.clstr | grep "\*" | awk -F"[>.]" '{print ">"$2}' | grep -Fwf - -A1 <(seqkit seq -w 0 SRR11015356_before_k75.unitigs.fa) | grep -v "^\-\-" > reps_unitigs_SRR11015356_before_k75.fa
```

### 1.5 Constructing cDBG @ k=75 for the representative sequences from step 1.4

```bash
bcalm -kmer-size 75 -max-memory 12000 -out reps_unitigs_SRR11015356_before_k75_after_k75.fa -in reps_unitigs_SRR11015356_before_k75.fa  -abundance-min 1
```

### 1.6 Extract the longest isoforms

```bash
python extract_longest_isoforms.py extractedGTF_gencode.v33.transcripts.fa
```

### 1.7 Remove overlapped transcripts using the GTF positions

```bash
python filter_overlaps.py longIso_extractedGTF_gencode.v33.transcripts.fa gencode.v33.basic.annotation.gtf
```

### 1.8 Filterout pseudogenes

```bash
cat nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa | seqkit grep -w 0 -v -n -r -p pseudogene > noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
```

> Reduced by 11205 sequences

```tsv
file                                                                 format  type  num_seqs      sum_len  min_len  avg_len  max_len
nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa           FASTA   DNA     46,934   97,352,548        8  2,074.2  205,012
noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa  FASTA   DNA     35,729   89,088,742        8  2,493.5  205,012
```


---

## 2. Alignment

### 2.1 Indexing

```bash
samtools faidx noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
bowtie2-build noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa noPsuedo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts
```

### 2.2 Before75_unitigs Alignment

```bash

# Alignment
bowtie2 -p 4 -x noPsuedo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts -f SRR11015356_before_k75.unitigs.fa -S bowtie2_SRR11015356_before_k75.unitigs.sam

# Converting SAM to BAM
samtools view -S -b bowtie2_SRR11015356_before_k75.unitigs.sam -o noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam

# Sorting the BAM
samtools sort noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam -o sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam

# Indexing the sorted BAM
samtools index sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam

```

Alignment summary

```txt
11824622 reads; of these:
  11824622 (100.00%) were unpaired; of these:
    10716275 (90.63%) aligned 0 times
    1010683 (8.55%) aligned exactly 1 time
    97664 (0.83%) aligned >1 times
9.37% overall alignment rate
```

### 2.3 Before75_after75 Alignment

```bash

# Alignment
bowtie2 -p 4 -x noPsuedo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts -f reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa -S reps_unitigs_SRR11015356_before_k75_after_k75.sam

## 2.2 Converting to sorted BAM
samtools view -S -b reps_unitigs_SRR11015356_before_k75_after_k75.sam -o noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools sort noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam -o sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam

```

Alignment summary

```txt
1993007 reads; of these:
  1993007 (100.00%) were unpaired; of these:
    1905087 (95.59%) aligned 0 times
    84044 (4.22%) aligned exactly 1 time
    3876 (0.19%) aligned >1 times
4.41% overall alignment rate
```
