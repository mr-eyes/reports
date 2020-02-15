<<NOTES
- https://www.biostars.org/p/262240/#262362

- Bowtie2, on the other hand, is designed to map reads continuously to the indexed reference.
  It just happens to be particularly good at dealing with multi-mapping reads, which is why people like using it to align to the transcriptome.

However, if you're interested in doing quantification with a known transcriptome in a generally well-annotated organism like mouse, 
I'd recommend doing quantification of the transcripts directly (e.g. using Salmon). 
If you then want to do DE, this can be followed with something like tximport to get results into your favorite DE tool (e.g. DESeq2, EdgeR, limma-voom).
NOTES

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

## ----------------------------------------------------

# 1. Aligning the "reps_unitigs_SRR11015356_before_k75.fa" on the "extractedGTF_gencode.v33.transcripts.fa"

## 1.1 Indexing 
#bowtie2-build /home/mabuelanin/Desktop/kexpression_experiment/symbolic/reports/report5_gencodeReprs/extractedGTF_gencode.v33.transcripts.fa bowtie2_extractedGTF_gencode.v33.transcripts
#samtools faidx /home/mabuelanin/Desktop/kexpression_experiment/symbolic/reports/report5_gencodeReprs/extractedGTF_gencode.v33.transcripts.fa # Important for IGV Visualization

## 1.1 Aligning
bowtie2 -x bowtie2_extractedGTF_gencode.v33.transcripts -f reps_unitigs_SRR11015356_before_k75.fa -S bowtie2_reps_unitigs_SRR11015356_before_k75.sam

## 1.2 Converting to sorted BAM
samtools view -S -b bowtie2_reps_unitigs_SRR11015356_before_k75.sam -o bowtie2_reps_unitigs_SRR11015356_before_k75.bam
samtools sort bowtie2_reps_unitigs_SRR11015356_before_k75.bam -o sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam

## 1.3 Generating coverage histogram
samtools index sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam
samtools depth sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam > sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam.cov
cat sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam.cov.hist

## 1.4 [Optional] View the alignment in the terminal
# samtools tview sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam extractedGTF_gencode.v33.transcripts.fa

## ---------------------------------------------------

# 2. Aligning the "reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa" on the "extractedGTF_gencode.v33.transcripts.fa"

## 2.1 Indexing [Done in step #1]
# bowtie2-build extractedGTF_gencode.v33.transcripts.fa bowtie2_extractedGTF_gencode.v33.transcripts
# samtools faidx extractedGTF_gencode.v33.transcripts.fa # Important for IGV Visualization


## 2.1 Aligning
bowtie2 -x bowtie2_extractedGTF_gencode.v33.transcripts -f reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa -S reps_unitigs_SRR11015356_before_k75_after_k75.sam


## 2.2 Converting to sorted BAM
samtools view -S -b reps_unitigs_SRR11015356_before_k75_after_k75.sam -o bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools sort bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam -o sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam


## 2.3 Generating coverage histogram
samtools index sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam
samtools depth sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam > sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov
cat sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov | awk -F'\t' '{print $3}' | sort -n | uniq -c > sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov.hist


## 2.4 [Optional] View the alignment in the terminal
# samtools tview sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam extractedGTF_gencode.v33.transcripts.fa


## ---------------------------------------------------

# 3. Merging the two histograms

paste sorted_bowtie2_reps_unitigs_SRR11015356_before_k75.bam.cov.hist sorted_bowtie2_reps_unitigs_SRR11015356_before_k75_after_k75.bam.cov.hist > before_after_cov.hist

# 4. IGV Visualization

## 4.1 Indexing the unitigs files


<<ALIGN_STATS
1.
2528169 reads; of these:
  2528169 (100.00%) were unpaired; of these:
    2381769 (94.21%) aligned 0 times
    38856 (1.54%) aligned exactly 1 time
    107544 (4.25%) aligned >1 times
5.79% overall alignment rate
----
2.
1993007 reads; of these:
  1993007 (100.00%) were unpaired; of these:
    1895443 (95.10%) aligned 0 times
    25712 (1.29%) aligned exactly 1 time
    71852 (3.61%) aligned >1 times
4.90% overall alignment rate
ALIGN_STATS