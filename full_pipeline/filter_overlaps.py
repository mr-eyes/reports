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