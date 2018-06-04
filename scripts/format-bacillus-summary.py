#! /usr/bin/env python3
"""
Read BCG and Bacillus (nonBCG) summary outputs and format them into a tabular
output for downstream analysis.

Example Summary Entry, 3 lines then empty line:
1. Bacillus thuringiensis YBT-1518, complete genome
6,002,284 bp circular DNA
NC_022873.1 GI:558678294

2. Bacillus cereus strain FORC60 chromosome, complete genome
5,361,178 bp circular DNA
NZ_CP020383.1 GI:1372048906

"""
import argparse as ap


def print_summary(entry, is_bcg):
    """Print parsed summary entry."""
    entry[0] = entry[0].replace("[", "").replace("]", "")
    species = entry[0].split(' ')[2]

    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(
        entry[2].split()[0].split('.')[0],  # accession
        entry[2].split()[1].replace('GI:', ''),  # gi
        is_bcg,
        True if species == 'anthracis' else False,  # is_ba
        species,
        entry[1].split(' ')[0],  # genome size
        entry[0].split(' ', 1)[1]  # description
    ))


def parse_summary(summary_file, is_bcg=False):
    """Parse summary outputs."""
    with open(summary_file, 'r') as summary_handle:
        entry = []
        for line in summary_handle:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            elif line:
                entry.append(line)
            else:
                # Parse the entry
                # entry[0] 1. Bacillus thuringiensis YBT-1518, complete genome
                # entry[1] 6,002,284 bp circular DNA
                # entry[2] NC_022873.1 GI:558678294
                print_summary(entry, is_bcg)
                entry.clear()
        print_summary(entry, is_bcg)


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='format-bacillus-summary.py', conflict_handler='resolve',
        description="Merge bcg and bacillus NCBI summaries into singel file."
    )

    parser.add_argument(
        'bcg', type=str, metavar="BCG_SUMMARY",
        help='NCBI summary output from BCG query.'
    )
    parser.add_argument(
        'bacillus', type=str, metavar="BACILLUS_SUMMARY",
        help='NCBI summary output from Bacillus (nonBCG) query.'
    )

    args = parser.parse_args()

    # Parse BCG summary
    print("#accession\tgi\tis_bcg\tis_ba\tspecies\tgenome_size\tdescription")
    parse_summary(args.bcg, is_bcg=True)
    parse_summary(args.bacillus)
