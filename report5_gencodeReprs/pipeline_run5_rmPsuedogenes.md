# Pipeline run3

## TODO

- [x] Remove the contained genes
- [ ] Remove pseudogenes
- [ ] Clustering CDHIT at 99% then repeat repor6.

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

## Filterout the overlapped shortest transcripts

```python
import sys
from itertools import groupby
import textwrap
import os


"""
This script take the transcripts records from GTF and the longest Isoforms fasta file which includes the longest isoform only from the same gene.
The output will be a new fasta file with the prefix "nonoverlapped_" for the longest nonoverlapped transcripts.
"""

def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

if len(sys.argv) < 3:
    sys.exit("run: python filter_overlaps.py <longIso_fasta> <transcripts_GTF>")

longIso_fasta = sys.argv[1]
transcripts_gtf = sys.argv[2]


longest_transcripts = set()
with open(longIso_fasta, 'r') as longIso_reader:
    for line in longIso_reader:
        if line[0] == ">":
            tr_id = line[1:].split('|')[0]
            longest_transcripts.add(tr_id)



interval_to_transcript = dict()
with open(transcripts_gtf, 'r') as gtf_reader:
    for line in gtf_reader:
        line = line.strip().split('\t')
        interval = (int(line[3]), int(line[4]))
        transcript_id = line[8].split(';')[1][16:-1]
        if transcript_id in longest_transcripts:
            interval_to_transcript[interval] = transcript_id


overlapped_transcripts = list()
filtered_transcripts = list()
prev_start, prev_end = next(iter(interval_to_transcript))
for interval, tr in interval_to_transcript.items():
    start = interval[0]
    end = interval[1]

    if ( start < prev_end ) or (prev_start > start):
        detected_overlap = [(prev_start,prev_end) , (start, end)]

        if (prev_end - prev_start) < (end - start):
            filtered_transcripts.append(interval_to_transcript[(prev_start, prev_end)])
        else:
            filtered_transcripts.append(interval_to_transcript[(start, end)])

        overlapped_transcripts.append(detected_overlap)

    prev_start, prev_end = start, end


assert len(filtered_transcripts) == len(overlapped_transcripts)

new_fasta = "nonoverlap_" + os.path.basename(longIso_fasta)
with open(new_fasta, 'w') as fasta_writer:
    for header, seq in fasta_iter(longIso_fasta):
        transcript_id = header.split('|')[0]
        if transcript_id not in filtered_transcripts:
            fasta_writer.write(f">{header}\n")
            fasta_writer.write(textwrap.fill(seq, 60) + '\n')
```


## **NEW** Filterout pseudogenes

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

## **1. Alignment**

```bash
# # 1. Aligning the "SRR11015356_before_k75.unitigs.fa" on the "noPseudo_longIso_extractedGTF_gencode.v33.transcripts.fa"

# ## 1.1 Indexing
bowtie2-build noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa noPsuedo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts
samtools faidx noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa  # Important for IGV Visualization

## 1.1 Aligning
bowtie2 -p 4 -x noPsuedo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts -f SRR11015356_before_k75.unitigs.fa -S bowtie2_SRR11015356_before_k75.unitigs.sam

## 1.2 Converting to sorted BAM
samtools view -S -b bowtie2_SRR11015356_before_k75.unitigs.sam -o noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam
samtools sort noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam -o sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam

## 1.3 Indexing BAM file
samtools index sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam

## 1.4 [Optional] View the alignment in the terminal
# samtools tview sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa

```

### 1.5 Alignment summary

```txt
11824622 reads; of these:
  11824622 (100.00%) were unpaired; of these:
    10716275 (90.63%) aligned 0 times
    1010683 (8.55%) aligned exactly 1 time
    97664 (0.83%) aligned >1 times
9.37% overall alignment rate
```

```bash

## ---------------------------------------------------

# 2. Aligning the "reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa" on the "nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa"

## 2.1 Indexing [Done in step #1]


## 2.1 Aligning
bowtie2 -p 4 -x noPsuedo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts -f reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa -S reps_unitigs_SRR11015356_before_k75_after_k75.sam


## 2.2 Converting to sorted BAM
samtools view -S -b reps_unitigs_SRR11015356_before_k75_after_k75.sam -o noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools sort noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam -o sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam


## 2.3 Generating coverage histogram
samtools index sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam


## 2.4 [Optional] View the alignment in the terminal
# samtools tview sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam extractedGTF_gencode.v33.transcripts.fa

```

## 2.5 Alginemnt summary

```txt
1993007 reads; of these:
  1993007 (100.00%) were unpaired; of these:
    1905087 (95.59%) aligned 0 times
    84044 (4.22%) aligned exactly 1 time
    3876 (0.19%) aligned >1 times
4.41% overall alignment rate
```

---

## **3. Get the primary alignments only from the BAM file**

> <https://www.biostars.org/p/259963/>

```bash

# Remove non-primary reads
samtools view -h -F 256 sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam > primary_sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam
samtools view -h -F 256 sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam > primary_sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam

# Convert SAM to BAM
# samtools view -S -b primary_bowtie2_SRR11015356_before_k75.unitigs.sam -o primary_bowtie2_SRR11015356_before_k75.unitigs.bam
# samtools view -S -b primary_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.sam -o primary_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam

# Sorting
samtools sort primary_sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam -o sorted_primary_noPsuedo_SRR11015356_before_k75.unitigs.bam
samtools sort primary_sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam -o sorted_primary_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam

# Indexing
samtools index sorted_primary_noPsuedo_SRR11015356_before_k75.unitigs.bam
samtools index sorted_primary_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam

```

## 4. Calculating depth and coverage

```bash

# Before 75
BAM_FILE1=sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam
BAM_FILE2=sorted_primary_noPsuedo_SRR11015356_before_k75.unitigs.bam
samtools depth ${BAM_FILE1} > ${BAM_FILE1}.cov
samtools depth ${BAM_FILE2} > ${BAM_FILE2}.cov
cat ${BAM_FILE1}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE1}.cov.hist
cat ${BAM_FILE2}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE2}.cov.hist

# After 75
BAM_FILE1=sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam
BAM_FILE2=sorted_primary_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools depth ${BAM_FILE1} > ${BAM_FILE1}.cov
samtools depth ${BAM_FILE2} > ${BAM_FILE2}.cov
cat ${BAM_FILE1}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE1}.cov.hist
cat ${BAM_FILE2}.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > ${BAM_FILE2}.cov.hist

# Merging histograms
paste sorted_primary_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov.hist  sorted_primary_noPsuedo_SRR11015356_before_k75.unitigs.bam.cov.hist > primary_before_after_cov.hist.tsv
paste sorted_noPseudo_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov.hist  sorted_noPsuedo_bowtie2_SRR11015356_before_k75.unitigs.bam.cov.hist > before_after_cov.hist.tsv

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

![](./run_5_IGV/1.png?raw=true)


## **6.1 Coverage**

**[Coverge](./run_5_coverage)**

## **6.2 IGV**

**[IGV](./run_5_IGV)**
