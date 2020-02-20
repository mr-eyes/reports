# Pipeline run3

## **Generate cDBGs**

```bash
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

start=`date +%s`


DATA_DIR=/home/mabuelanin/Desktop/kexpression_experiment/symbolic/data/drtamer_data

echo "Generating cDBG for raw reads (two files) with k=75"
/usr/bin/time -v bcalm -kmer-size 75 -max-memory 12000 -out SRR11015356_before_k75 -in list_reads &> bcalm_before_75.log

echo "CD-HIT-EST Clustering for SRR11015356_before_k75"
/usr/bin/time -v cd-hit-est -i SRR11015356_before_k75.unitigs.fa -n 11 -c 0.95 -o clusters_SRR11015356_before_k75 -d 0 -T 0 -M 12000 &> cdhit_SRR11015356_before_k75.log

echo "Exporting representative sequences only from the CDHIT*clst and for SRR11015356_before_k75.unitigs.fa ..."
cat clusters_SRR11015356_before_k75.clstr | grep "\*" | awk -F"[>.]" '{print ">"$2}' | grep -Fwf - -A1 <(seqkit seq -w 0 SRR11015356_before_k75.unitigs.fa) | grep -v "^\-\-" > reps_unitigs_SRR11015356_before_k75.fa

echo "Creating cDBG for reps_unitigs_SRR11015356_before_k75.fa with k=25"
/usr/bin/time -v bcalm -kmer-size 25 -max-memory 12000 -out reps_unitigs_SRR11015356_before_k75_after_k25.fa -in reps_unitigs_SRR11015356_before_k75.fa  -abundance-min 1 &> bcalm_before_75_after_25.log

echo "Creating cDBG for reps_unitigs_SRR11015356_before_k75.fa with k=75"
/usr/bin/time -v bcalm -kmer-size 75 -max-memory 12000 -out reps_unitigs_SRR11015356_before_k75_after_k75.fa -in reps_unitigs_SRR11015356_before_k75.fa  -abundance-min 1 &> bcalm_before_75_after_75.log


end=`date +%s`

runtime=$((end-start))

echo "DONE SUCCESSFULLY in ${runtime}"

```

---

## **Extract the longest isoforms** > *longIso_extractedGTF_gencode.v33.transcripts.fa*

```python

import sys
from itertools import groupby
import textwrap
import os

"""
Extracting the longest isoforms from the CHR Basic annotation from gencode and the whole transcriptome file from Gencode.
"""

def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


fasta_file = "extractedGTF_gencode.v33.transcripts.fa"
gene_to_transcripts = dict()

for header, seq in fasta_iter(fasta_file):
    header_info = header.split('|')
    transcript_id = header_info[0]
    gene_id = header_info[1]
    seq_len = len(seq)
    if gene_id not in gene_to_transcripts:
        gene_to_transcripts[gene_id] = dict()

    gene_to_transcripts[gene_id][transcript_id] = seq_len

longest_isoforms = set()

for gene, transcripts in gene_to_transcripts.items():
    longest_isoforms.add(max(transcripts, key=transcripts.get))

new_fasta = "longIso_" + os.path.basename(fasta_file)

with open(new_fasta, 'w') as fasta_writer:
    for header, seq in fasta_iter(fasta_file):
        transcript_id = header.split('|')[0]
        if transcript_id in longest_isoforms:
            fasta_writer.write(f">{header}\n")
            fasta_writer.write(textwrap.fill(seq, 60) + '\n')

```

---

## **1. Alignment**

```bash
# # 1. Aligning the "SRR11015356_before_k75.unitigs.fa" on the "longIso_extractedGTF_gencode.v33.transcripts.fa"

# ## 1.1 Indexing
bowtie2-build /home/mabuelanin/Desktop/kexpression_experiment/symbolic/reports/report5_gencodeReprs/longIso_extractedGTF_gencode.v33.transcripts.fa longIso_extractedGTF_gencode.v33.transcripts
samtools faidx /home/mabuelanin/Desktop/kexpression_experiment/symbolic/reports/report5_gencodeReprs/longIso_extractedGTF_gencode.v33.transcripts.fa  # Important for IGV Visualization

## 1.1 Aligning
bowtie2 -x longIso_extractedGTF_gencode.v33.transcripts -f SRR11015356_before_k75.unitigs.fa -S bowtie2_SRR11015356_before_k75.unitigs.sam

## 1.2 Converting to sorted BAM
samtools view -S -b bowtie2_SRR11015356_before_k75.unitigs.sam -o bowtie2_SRR11015356_before_k75.unitigs.bam
samtools sort bowtie2_SRR11015356_before_k75.unitigs.bam -o sorted_bowtie2_SRR11015356_before_k75.unitigs.bam

## 1.3 Indexing BAM file
samtools index sorted_bowtie2_SRR11015356_before_k75.unitigs.bam

## 1.4 [Optional] View the alignment in the terminal
# samtools tview sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam longIso_extractedGTF_gencode.v33.transcripts.fa

```

### 1.5 Alignment summary

```txt
11824622 reads; of these:
  11824622 (100.00%) were unpaired; of these:
    10621664 (89.83%) aligned 0 times
    929026 (7.86%) aligned exactly 1 time
    273932 (2.32%) aligned >1 times
10.17% overall alignment rate
```

```bash

## ---------------------------------------------------

# 2. Aligning the "reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa" on the "longIso_extractedGTF_gencode.v33.transcripts.fa"

## 2.1 Indexing [Done in step #1]


## 2.1 Aligning
bowtie2 -x longIso_extractedGTF_gencode.v33.transcripts -f reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa -S reps_unitigs_SRR11015356_before_k75_after_k75.sam


## 2.2 Converting to sorted BAM
samtools view -S -b reps_unitigs_SRR11015356_before_k75_after_k75.sam -o bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools sort bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam -o sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam


## 2.3 Generating coverage histogram
samtools index sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam


## 2.4 [Optional] View the alignment in the terminal
# samtools tview sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam extractedGTF_gencode.v33.transcripts.fa

```

## 2.5 Alginemnt summary

```txt
1993007 reads; of these:
  1993007 (100.00%) were unpaired; of these:
    1899598 (95.31%) aligned 0 times
    80909 (4.06%) aligned exactly 1 time
    12500 (0.63%) aligned >1 times
4.69% overall alignment rate
```

---


## **3. Get the primary alignments only from the BAM file**

> <https://www.biostars.org/p/259963/>

```bash

samtools view -F 256 sorted_bowtie2_SRR11015356_before_k75.unitigs.bam > primary_sorted_bowtie2_SRR11015356_before_k75.unitigs.bam
samtools view -F 256 sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam > primary_sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools index primary_sorted_bowtie2_SRR11015356_before_k75.unitigs.bam
samtools index primary_sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam

```

## 4. Calculating depth and coverage

```bash

# Before 75

BAM_FILE1=sorted_bowtie2_SRR11015356_before_k75.unitigs.bam
BAM_FILE2=primary_sorted_bowtie2_SRR11015356_before_k75.unitigs.bam
samtools depth ${BAM_FILE1} > ${BAM_FILE1}.cov
samtools depth ${BAM_FILE2} > ${BAM_FILE2}.cov
cat ${BAM_FILE1}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE1}.cov.hist
cat ${BAM_FILE2}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE2}.cov.hist

# After 75
BAM_FILE1=primary_sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
BAM_FILE2=sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools depth ${BAM_FILE1} > ${BAM_FILE1}.cov
samtools depth ${BAM_FILE2} > ${BAM_FILE2}.cov
cat ${BAM_FILE1}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE1}.cov.hist
cat ${BAM_FILE2}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE2}.cov.hist

# Merging histograms
paste primary_sorted_bowtie2_SRR11015356_before_k75.unitigs.bam.cov.hist primary_sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov.hist > primary_before_after_cov.hist.tsv
paste sorted_bowtie2_SRR11015356_before_k75.unitigs.bam.cov.hist sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov.hist > before_after_cov.hist.tsv

```

## 5. Delete dead files

```bash

rm bowtie2_SRR11015356_before_k75.unitigs.bam
rm bowtie2_SRR11015356_before_k75.unitigs.sam
rm bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
rm bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.sam

```

--

# **6. Results**

## **6.1 Coverage**

**[Coverge](./run_3_coverage)**

## **6.2 IGV**

**[IGV](./run_3_IGV)**
