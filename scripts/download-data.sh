#! /bin/bash
# download-sequences.sh
# This script automates the download of completed genomes and rRNAs.
TOP_DIR=$(pwd)
SEQ_DIR=$1

mkdir -p ${SEQ_DIR}
# Completed Genomes
mkdir ${SEQ_DIR}/completed-genomes
entrez-download.py /opt/data/bacillus-genomes.txt ${SEQ_DIR}/completed-genomes

# Bacillus rRNA
entrez-download.py /opt/data/bacillus-rrna.txt ${SEQ_DIR}/ncbi-rrna.fasta --multi_fasta

# Silva rRNA
SILVA_URL='https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports'
LSUParc="SILVA_132_LSUParc_tax_silva.fasta.gz"
SSUParc="SILVA_132_SSUParc_tax_silva.fasta.gz"
if [ ! -f ${SEQ_DIR}/${LSUParc} ]; then
    wget -O ${SEQ_DIR}/${LSUParc} ${SILVA_URL}/${LSUParc}
else
    echo "Skipping ${SEQ_DIR}/${LSUParc}, already exists..."
fi

if [ ! -f ${SEQ_DIR}/${SSUParc} ]; then
    wget -O ${SEQ_DIR}/${SSUParc} ${SILVA_URL}/${SSUParc}
else
    echo "Skipping ${SEQ_DIR}/${SSUParc}, already exists..."
fi

if [ ! -f ${SEQ_DIR}/rrna.fasta ]; then
    cp ${SEQ_DIR}/ncbi-rrna.fasta ${SEQ_DIR}/rrna.fasta
    zcat ${SEQ_DIR}/${SSUParc} ${SEQ_DIR}/${SSUParc} >> ${SEQ_DIR}/rrna.fasta
else
    echo "Skipping ${SEQ_DIR}/rrna.fasta, already exists..."
fi


# Count 31mers for each FASTA file and gzip the file
for fasta in `find ${SEQ_DIR} -name "*.fasta"`; do
    if [ ! -f ${fasta}.jf ]; then
        jellyfish count -C -m 31 -s 1M -o ${fasta}.jf ${fasta}
    else
        echo "Skipping ${fasta}.jf, already exists..."
    fi
done
