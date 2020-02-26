# Report 10 (Different aligners)

## Summary

- [x] 1. Try [LAST](http://last.cbrc.jp/doc/last-tutorial.html)

## **1 LAST**

```bash
conda install -c bioconda last
```

### 1.1 SRR11015356_before_k75.unitigs.fa

#### 1.1.2 Indexing

```bash
lastdb -v -P 4 -cR01 LAST_transcriptome noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
```

#### 1.1.3 Alignment

```bash
lastal -P 4 -v LAST_transcriptome SRR11015356_before_k75.unitigs.fa > LAST_SRR11015356_before_k75.maf
```

#### 1.1.3 Converting MAF to SAM [reference](http://last.cbrc.jp/doc/maf-convert.html)

```bash
maf-convert sam LAST_SRR11015356_before_k75.maf > LAST_SRR11015356_before_k75.sam
samtools view -bt noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa.fai LAST_SRR11015356_before_k75.sam > LAST_SRR11015356_before_k75.bam
samtools sort LAST_SRR11015356_before_k75.bam -o sorted_LAST_SRR11015356_before_k75.bam
samtools index sorted_LAST_SRR11015356_before_k75.bam
```

---

### 2.1 reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa

#### 2.1.2 ~Indexing~ Done

```bash
# lastdb -v -P 4 -cR01 LAST_transcriptome noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa
```

#### 2.1.3 Alignment

```bash
lastal -P 4 -v LAST_transcriptome reps_unitigs_SRR11015356_before_k75_after_k75.fa.unitigs.fa > LAST_reps_SRR11015356_before_k75_after_k75.maf
```

#### 2.1.3 Converting MAF to SAM

```bash
maf-convert sam LAST_reps_SRR11015356_before_k75_after_k75.maf > LAST_reps_SRR11015356_before_k75_after_k75.sam
samtools view -bt noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa.fai LAST_reps_SRR11015356_before_k75_after_k75.sam > LAST_reps_SRR11015356_before_k75_after_k75.bam
samtools sort LAST_reps_SRR11015356_before_k75_after_k75.bam -o sorted_LAST_reps_SRR11015356_before_k75_after_k75.bam
samtools index sorted_LAST_reps_SRR11015356_before_k75_after_k75.bam
```
