#! /usr/bin/env python3
"""Download FASTA files for a given set of NCBI accessions."""
import argparse as ap
import os
import sys
import time
from Bio import Entrez
Entrez.email = 'robert.petit@emory.edu'
DATABASE = 'nucleotide'
RETTYPE = 'fasta'
RETMODE = 'text'


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='entrez-download.py', conflict_handler='resolve',
        description="Download FASTA files for a given set of NCBI accessions."
    )

    parser.add_argument('accessions', type=str, metavar="QUERY",
                        help='File with NCBI accessions in the first column.')
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('Directory to write sequences to.'))
    parser.add_argument('--multi_fasta', action="store_true",
                        help='Write output as single multi-FASTA file.')

    args = parser.parse_args()
    output = args.output
    write_mode = 'w'
    if args.multi_fasta:
        write_mode = 'a'
        if os.path.exists(output):
            print('{0} exists, skipping downloads...'.format(output))
            sys.exit()

    # Parse BCG summary
    with open(args.accessions, 'r') as accession_handle:
        for line in accession_handle:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            else:
                fasta_output = output
                accession = line.split('\t')[0]
                if not args.multi_fasta:
                    fasta_output = '{0}/{1}.fasta'.format(output, accession)

                if os.path.exists(output) and not args.multi_fasta:
                    print("Skipping {0}".format(accession))
                else:
                    print("Downloading {0}".format(accession))
                    with open(fasta_output, write_mode) as fasta_handle:
                        efetch = Entrez.efetch(
                            db=DATABASE, id=accession, rettype=RETTYPE,
                            retmode=RETMODE
                        )
                        fasta_handle.write(efetch.read())
                        efetch.close()
                    time.sleep(0.5)
