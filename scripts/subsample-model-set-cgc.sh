#! /bin/bash
set -x
CPU=$2
BA_KMERS="/opt/results/ba-specific-kmers.fasta"
BCG_KMERS="/opt/results/bcg-specific-kmers.fasta"
LEF_KMERS="/opt/results/lef-specific-kmers.fasta"
SUBSAMPLE_COVERAGES="/opt/data/subsample-coverages.txt"

# B. anthracis sequencing
fastq-dump --split-spot -I --fasta 0 $1

# Remove sequences with Ns
grep -v "^>" $1.fasta > $1-sequence.txt

# Check group specific kmers on whole set
jellyfish count -C  -m 31 -s 5M -t 22 -o $1.jf $1.fasta
jellyfish query -s ${BA_KMERS} -o $1-ba.txt $1.jf
jellyfish query -s ${BCG_KMERS} -o $1-bcg.txt $1.jf
jellyfish query -s ${LEF_KMERS} -o $1-lef.txt $1.jf

# Determine lethal factor detection limit
mkdir temp_folder/
subsample-model-set.py $1 $1-sequence.txt ${SUBSAMPLE_COVERAGES} temp_folder/ \
    subsample/ ${BA_KMERS} ${BCG_KMERS} ${LEF_KMERS} --cpu ${CPU} --replicates 20

mkdir full-counts/
mv $1-ba.txt $1-bcg.txt $1-lef.txt full-counts/

# Summarize the Counts
find -mindepth 3 -maxdepth 3 -type d -path "*subsample*" | \
xargs -I {} -P ${CPU} -n 1 summarize-model-set.py {} summary/

# Tarball the results
tar -cf $1.tar subsample/ full-counts/ summary/
rm -rf $1.fasta  $1.jf  $1-sequence.txt temp_folder/ subsample/ full-counts/ summary/
