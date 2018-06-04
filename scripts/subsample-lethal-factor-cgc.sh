#! /bin/bash
set -x
BA_KMERS="/opt/results/ba-specific-kmers.fasta"
BCG_KMERS="/opt/results/bcg-specific-kmers.fasta"
LEF_KMERS="/opt/results/lef-specific-kmers.fasta"

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
subsample-lethal-factor.py $1 $1-sequence.txt temp_folder/ subsample-lethal/ \
    ${BA_KMERS} ${BCG_KMERS} ${LEF_KMERS} --cpu 30 --replicates 100

mkdir full-counts/
mv $1-ba.txt $1-bcg.txt $1-lef.txt full-counts/

# Tarball the results
tar -cf $1.tar subsample-lethal/ full-counts/
rm -rf $1.fasta  $1.jf  $1-sequence.txt temp_folder/ subsample-lethal/ full-counts/
