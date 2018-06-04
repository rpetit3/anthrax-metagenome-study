#! /usr/bin/env python3
"""Download FASTA files for a given set of NCBI accessions."""
import argparse as ap
from collections import Counter, OrderedDict
import os
from itertools import islice
import sys
import time
from Bio import Entrez
Entrez.email = 'robert.petit@emory.edu'
DATABASE = 'nucleotide'
RETTYPE = 'fasta_cds_na'
RETMODE = 'text'


def download_plasmids(accessions, output):
    sequences = OrderedDict()
    with open(accessions, 'r') as accession_handle:
        for line in accession_handle:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            else:
                accession = line.split('\t')[0]
                fasta_output = '{0}/{1}.fasta'.format(output, accession)
                if os.path.exists(fasta_output):
                    print("Skipping {0}".format(accession))
                else:
                    print("Downloading {0}".format(accession))
                    with open(fasta_output, 'w') as fasta_handle:
                        efetch = Entrez.efetch(
                            db=DATABASE, id=accession, rettype=RETTYPE,
                            retmode=RETMODE
                        )
                        fasta_handle.write(efetch.read())
                        efetch.close()
                    time.sleep(0.5)
                sequences[accession] = read_fasta(fasta_output)
    return sequences


def read_fasta(fasta):
    """Return a list of seqeunces from a given FASTA file."""
    if os.path.exists(fasta):
        header = None
        seq = []
        records = {}
        with open(fasta, 'r') as fasta_handle:
            for line in fasta_handle:
                line = line.rstrip()
                if line.startswith('>'):
                    if seq:
                        records[header] = ''.join(seq)
                        seq = []
                    header = line
                else:
                    seq.append(line)

        records[header] = ''.join(seq)

        return records


def split_into_kmers(seq, k=31):
    """
    Source: http://stackoverflow.com/questions/7636004/                   \
                 python-split-string-in-moving-window
    Returns a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn),
    """
    iterator = iter(seq)
    result = tuple(islice(iterator, k))
    if len(result) == k:
        yield result
    for elem in iterator:
        result = result[1:] + (elem,)
        yield result


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='process-lethal-factor.py', conflict_handler='resolve',
        description="Pull out lethal factors from pXO1 plasmids."
    )

    parser.add_argument('accessions', type=str, metavar="QUERY",
                        help='File with pXO1 NCBI accessions.')
    parser.add_argument('genomes', type=str, metavar="OUTPUT",
                        help=('Directory to write sequences to.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('Directory to final kmer sequences to.'))

    args = parser.parse_args()
    plasmids = download_plasmids(args.accessions, args.genomes)
    common_kmers = Counter()
    for plasmid, records in plasmids.items():
        for record, sequence in records.items():
            if "lethal factor" in record:
                for i in split_into_kmers(sequence):
                    common_kmers["".join(i)] += 1

    final_output = "{0}/lef-specific-kmers.fasta".format(args.output)
    with open(final_output, 'w') as fasta:
        for kmer, count in common_kmers.items():
            fasta.write(">{0}-{1}\n{0}\n".format(kmer, count))
