# **Report 6: The Optimum K-mer**

## *Summary*

1. Abundance plot:
    1. Reference: `longIso_extractedGTF_gencode.v33.transcripts.fa`
    2. k-mers to test: {25,31,41,51,75,81}
    3. Scatter plot with x-axis:kmer_size, y-axis:number of unique kmers with abundance > 1.
    4. By visual inspection, select the best k-mer size for the next steps.
2. Create cDBG with the optimum kmer size from step #1.
3. Clustering similarity threshold:
    1. Run the CD-HIT-EST on the unitigs at different similarity threshold {90,93,96,99}.
    2. Scatter plot with x-axis: similarity threshold, y-axis: number of clusters with size > 1.
    