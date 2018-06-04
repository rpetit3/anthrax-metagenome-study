#! /usr/bin/env python3
"""Create random subsampels of input sequences."""
import argparse as ap
import gzip
import glob
from multiprocessing import Pool
import os
import random
import subprocess
import numpy as np
GENOME_SIZE = None
BA_LENGTH = None
BA_SEQUENCES = None
BA_TOTAL = None
BCG_LENGTH = None
BCG_SEQUENCES = None
BCG_TOTAL = None
BA_KMERS = None
BCG_KMERS = None
LEF_KMERS = None


def read_coverages(input_file):
    """Return coverage values."""
    coverages = []
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            coverages.append(float(line.rstrip()))
    return coverages


def read_sequences(input_file, min_length=None):
    """Return lines in a text file as a list."""
    lines = []
    lengths = []
    total = 1
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            line = line.rstrip()
            length = len(line)
            if min_length:
                if length >= min_length:
                    lines.append(line)
                    if total <= 100000:
                        lengths.append(length)
                    total += 1
            else:
                lines.append(line)
                lengths.append(length)
                if total <= 100000:
                    lengths.append(length)

    length_stats = get_coverage_stats(lengths)
    return [lines, int(length_stats['mean'])]


def get_coverage_stats(coverage):
    """Return summary stats of a set of coverages."""
    np_array = np.array(coverage)
    return {
        'min': min(coverage) if coverage else 0,
        'median': int(np.median(np_array)) if coverage else 0,
        'mean': np.mean(np_array) if coverage else 0,
        'max': max(coverage) if coverage else 0
    }


def output_handler(output, redirect='>'):
    if output:
        return [open(output, 'w'), '{0} {1}'.format(redirect, output)]
    else:
        return [subprocess.PIPE, '']


def run_command(cmd, cwd=os.getcwd(), stdout=False, stderr=False, shell=False):
    """Execute a single command and return STDOUT and STDERR."""
    stdout, stdout_str = output_handler(stdout)
    stderr, stderr_str = output_handler(stderr, redirect='2>')

    p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, cwd=cwd,
                         shell=shell)

    return p.communicate()


def subsample(opts):
    """Subsample coverages."""
    prefix = opts[0]
    bcg_coverage = opts[1]
    ba_coverage = opts[2]
    replicate = opts[3]
    working_dir = './'
    output_dir = './'
    bcg_reads = 0
    ba_reads = 0
    basename = "replicate-{0}-{1:03d}".format(prefix, replicate)

    if bcg_coverage or ba_coverage:
        if not os.path.exists('{0}/{1}-lef.txt.gz'.format(output_dir, basename)):
            fasta = []
            fasta_output = "{0}/{1}.fasta".format(working_dir, basename)
            random_seed = None

            if bcg_coverage and ba_coverage:
                bcg_reads = int(GENOME_SIZE * float(bcg_coverage) / BCG_LENGTH)
                ba_reads = int(GENOME_SIZE * float(ba_coverage) / BA_LENGTH)
                random_seed = (
                    int(bcg_coverage * 100) * int(ba_coverage * 1000) * replicate + bcg_reads + ba_reads
                )
            elif bcg_coverage:
                bcg_reads = int(GENOME_SIZE * float(bcg_coverage) / BCG_LENGTH)
                random_seed = (
                    int(bcg_coverage * 100) * replicate + bcg_reads
                )
            else:
                ba_reads = int(GENOME_SIZE * float(ba_coverage) / BA_LENGTH)
                random_seed = (
                    int(ba_coverage * 1000) * replicate + ba_reads
                )

            if bcg_coverage:
                bcg_reads = int(GENOME_SIZE * float(bcg_coverage) / BCG_LENGTH)
                random.seed(random_seed)
                for element in random.sample(range(BCG_TOTAL), bcg_reads):
                    fasta.append(">{0}\n".format(element))
                    fasta.append("{0}\n".format(BCG_SEQUENCES[element]))

            if ba_coverage:
                ba_reads = int(GENOME_SIZE * float(ba_coverage) / BA_LENGTH)
                random.seed(random_seed)
                for element in random.sample(range(BA_TOTAL), ba_reads):
                    fasta.append(">{0}\n".format(element))
                    fasta.append("{0}\n".format(BA_SEQUENCES[element]))

            print("\tReplicate: {0} Random Seed: {1} Reads: BCG {2} BA {3}".format(
                replicate, random_seed, bcg_reads, ba_reads
            ))

            with open(fasta_output, 'w') as fasta_handle:
                fasta_handle.write("".join(fasta))

            # Count kmers
            jellyfish = '{0}/{1}.jf'.format(working_dir, basename)
            run_command(['jellyfish', 'count', '-C', '-t', '4', '-m', '31',
                         '-s', '5M', '-o', jellyfish, fasta_output])
            run_command(['rm', fasta_output])
            ba_txt = '{0}/{1}-ba.txt'.format(working_dir, basename)
            run_command(
                ['jellyfish', 'query', '-s', BA_KMERS, '-o', ba_txt, jellyfish]
            )
            run_command(['gzip', '-f', ba_txt])
            run_command(['mv', '{0}.gz'.format(ba_txt), output_dir])

            bcg_txt = '{0}/{1}-bcg.txt'.format(working_dir, basename)
            run_command(
                ['jellyfish', 'query', '-s', BCG_KMERS, '-o', bcg_txt, jellyfish]
            )
            run_command(['gzip', '-f', bcg_txt])
            run_command(['mv', '{0}.gz'.format(bcg_txt), output_dir])

            lef_txt = '{0}/{1}-lef.txt'.format(working_dir, basename)
            run_command(
                ['jellyfish', 'query', '-s', LEF_KMERS, '-o', lef_txt, jellyfish]
            )
            run_command(['gzip', '-f', lef_txt])
            run_command(['mv', '{0}.gz'.format(lef_txt), output_dir])
            run_command(['rm', jellyfish])


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='subsample-sequences-cgc.py', conflict_handler='resolve',
        description="Create random subsamples of input sequences."
    )

    parser.add_argument('ba_sequences', type=str, metavar="BA_SEQS",
                        help='File of B. anthracis sequences, one per line.')
    parser.add_argument('bcg_sequences', type=str, metavar="BCG_SEQS",
                        help='File of B. cereus sequences, one per line.')
    parser.add_argument('ba', type=str, metavar="BA_KMERS",
                        help=('BA specific kmers.'))
    parser.add_argument('bcg', type=str, metavar="BCG_KMERS",
                        help=('BCG specific kmers.'))
    parser.add_argument('lef', type=str, metavar="LEF_KMERS",
                        help=('Lethal factor kmers.'))
    parser.add_argument('ba_coverage', type=float, metavar="BA_COVERAGE",
                        help=('Coverage to subsample B. anthracis to.'))
    parser.add_argument('bcg_coverage', type=float, metavar="BCG_COVERAGE",
                        help=('Coverage to subsample B. cereus to.'))
    parser.add_argument('--genome_size', metavar="INT", type=int,
                        default=5200000,
                        help='Genome size (Default 5.2Mb)')
    parser.add_argument('--length', metavar="INT", type=int, default=100,
                        help='Per line sequence length (Default 100)')
    parser.add_argument('--replicates', metavar="INT", type=int,
                        default=20,
                        help='Number of replicates per coverage (Default 20)')
    parser.add_argument('--cpu', metavar="INT", type=int, default=1,
                        help='Total number of processes to launch (Default 1)')

    args = parser.parse_args()

    BA_SEQUENCES, BA_LENGTH = read_sequences(args.ba_sequences,
                                             min_length=args.length)
    BCG_SEQUENCES, BCG_LENGTH = read_sequences(args.bcg_sequences,
                                               min_length=args.length)

    print("Mean Read Lengths: BCG {0}bp, BA {1}bp".format(
        BCG_LENGTH, BA_LENGTH
    ))

    BA_TOTAL = len(BA_SEQUENCES)
    BCG_TOTAL = len(BCG_SEQUENCES)
    BA_KMERS = args.ba
    BCG_KMERS = args.bcg
    LEF_KMERS = args.lef
    GENOME_SIZE = args.genome_size

    bcg_coverage = float(args.bcg_coverage)
    ba_coverage = float(args.ba_coverage)
    basename = "{0}-{1}".format(bcg_coverage, ba_coverage)
    for replicate in range(1, args.replicates + 1):
        subsample([basename, bcg_coverage, ba_coverage, replicate])

    tar_cmd = ['tar', '-cvf', 'subsample-{0}.tar'.format(basename)]
    rm_cmd = ['rm']
    for count in glob.glob("./*.txt.gz"):
        tar_cmd.append(count)
        rm_cmd.append(count)
    run_command(tar_cmd)
    run_command(rm_cmd)
