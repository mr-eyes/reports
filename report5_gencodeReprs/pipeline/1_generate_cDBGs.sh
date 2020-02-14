# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

start=`date +%s`


DATA_DIR=/home/mabuelanin/Desktop/kexpression_experiment/symbolic/data/drtamer_data

echo "Generating cDBG for raw reads (two files) with k=25"
/usr/bin/time -v bcalm -kmer-size 25 -max-memory 12000 -out SRR11015356_before_k25 -in list_reads &> bcalm_before_25.log

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
