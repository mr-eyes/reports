import kProcessor as kp

fasta_file="/home/mabuelanin/Desktop/kexpression_experiment/symbolic/data/drtamer_data/SRR11015356_k75.unitigs.fa"
names_file=fasta_file + ".names"


kSize = 25
Q = 33
mode = 1 # kmers
chunk_size = int(1e3)

#kf = kp.kDataFrameMQF(kSize, Q, mode)

# kf = kp.kDataFrameBMQF(kSize)

kf = kp.kDataFramePHMAP(kSize)

print("Indexing ...")

ckf = kp.index(kf, {"mode" : mode}, fasta_file , chunk_size, names_file)

print("Serializing the index ...")

ckf.save("idx_k25_SRR11015356_k75.unitigs")
