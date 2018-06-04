#! /usr/bin/env python3
"""Create random subsampels of input sequences."""
import argparse as ap
import glob
import gzip
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


def has_hit(counts):
    with open(counts, 'r') as count_handle:
        for line in count_handle:
            if int(line.rstrip().split()[1]):
                return True
    return False


def subsample(opts):
    """Subsample coverages."""
    working_dir = opts[0]
    output_dir = opts[1]
    bcg_coverage = opts[2]
    ba_coverage = opts[3]
    replicate = opts[4] + 100
    bcg_reads = 0
    ba_reads = 0
    basename = "replicate-{0:03d}".format(replicate)

    if bcg_coverage or ba_coverage:
        if not os.path.exists('{0}/{1}-lef.txt.gz'.format(output_dir, basename)):
            fasta = []
            fasta_output = "{0}/{1}.fasta".format(working_dir, basename)
            random_seed = None

            if bcg_coverage and ba_coverage:
                bcg_reads = int(GENOME_SIZE * float(bcg_coverage) / BCG_LENGTH)
                ba_reads = int(GENOME_SIZE * float(ba_coverage) / BA_LENGTH)
                random_seed = (
                    int(bcg_coverage * 100) * int(ba_coverage * 10000) * replicate + bcg_reads + ba_reads
                )
            elif bcg_coverage:
                bcg_reads = int(GENOME_SIZE * float(bcg_coverage) / BCG_LENGTH)
                random_seed = (
                    int(bcg_coverage * 100) * replicate + bcg_reads
                )
            else:
                ba_reads = int(GENOME_SIZE * float(ba_coverage) / BA_LENGTH)
                random_seed = (
                    int(ba_coverage * 10000) * replicate + ba_reads
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

            with open(fasta_output, 'w') as fasta_handle:
                fasta_handle.write("".join(fasta))
            run_command(['mv', fasta_output, output_dir])
            # Count kmers
            """
            jellyfish = '{0}/{1}.jf'.format(working_dir, basename)
            run_command(['jellyfish', 'count', '-C', '-m', '31',
                         '-s', '5M', '-o', jellyfish, fasta_output])
            run_command(['rm', fasta_output])
            ba_txt = '{0}/{1}-ba.txt'.format(working_dir, basename)
            run_command(
                ['jellyfish', 'query', '-s', BA_KMERS, '-o', ba_txt, jellyfish]
            )
            ba_hit = has_hit(ba_txt)

            bcg_txt = '{0}/{1}-bcg.txt'.format(working_dir, basename)
            run_command(
                ['jellyfish', 'query', '-s', BCG_KMERS, '-o', bcg_txt, jellyfish]
            )
            bcg_hit = has_hit(bcg_txt)

            lef_txt = '{0}/{1}-lef.txt'.format(working_dir, basename)
            run_command(
                ['jellyfish', 'query', '-s', LEF_KMERS, '-o', lef_txt, jellyfish]
            )
            run_command(['rm', jellyfish])

            if ba_hit and bcg_hit:
                print("\tSUCCESS: Replicate: {0} Random Seed: {1} Reads: BCG {2} BA {3}".format(
                    replicate, random_seed, bcg_reads, ba_reads
                ))
                run_command(['gzip', '-f', bcg_txt])
                run_command(['gzip', '-f', lef_txt])
                run_command(['gzip', '-f', ba_txt])
                run_command(['mv', '{0}.gz'.format(ba_txt), output_dir])
                run_command(['mv', '{0}.gz'.format(bcg_txt), output_dir])
                run_command(['mv', '{0}.gz'.format(lef_txt), output_dir])
            else:
                run_command(['rm', bcg_txt])
                run_command(['rm', lef_txt])
                run_command(['rm', ba_txt])
            """
        else:
            print("\tSkipping replicate: {0}, already completed".format(replicate))


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='subsample-sequences.py', conflict_handler='resolve',
        description="Create random subsamples of input sequences."
    )

    parser.add_argument('ba_sequences', type=str, metavar="BA_SEQS",
                        help='File of B. anthracis sequences, one per line.')
    parser.add_argument('bcg_sequences', type=str, metavar="BCG_SEQS",
                        help='File of B. cereus sequences, one per line.')
    parser.add_argument('ba_coverages', type=str, metavar="BA_COVERAGES",
                        help=('Coverages to subsample B. anthracis to.'))
    parser.add_argument('bcg_coverages', type=str, metavar="BCG_COVERAGES",
                        help=('Coverages to subsample B. cereus to.'))
    parser.add_argument('working', type=str, metavar="WORKING_DIR",
                        help=('Directory to put temporary files.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('Directory to subsampled FASTA files to.'))
    parser.add_argument('ba', type=str, metavar="BA_KMERS",
                        help=('BA specific kmers.'))
    parser.add_argument('bcg', type=str, metavar="BCG_KMERS",
                        help=('BCG specific kmers.'))
    parser.add_argument('lef', type=str, metavar="LEF_KMERS",
                        help=('Lethal factor kmers.'))
    parser.add_argument('--genome_size', metavar="INT", type=int,
                        default=5200000,
                        help='Genome size (Default 5.2Mb)')
    parser.add_argument('--length', metavar="INT", type=int, default=100,
                        help='Per line sequence length (Default 100)')
    parser.add_argument('--replicates', metavar="INT", type=int,
                        default=20,
                        help='Number of replicates per coverage (Default 100)')
    parser.add_argument('--cpu', metavar="INT", type=int, default=23,
                        help='Total number of processes to launch (Default 1)')

    args = parser.parse_args()

    BA_SEQUENCES, BA_LENGTH = read_sequences(args.ba_sequences,
                                             min_length=args.length)
    BCG_SEQUENCES, BCG_LENGTH = read_sequences(args.bcg_sequences,
                                               min_length=args.length)

    print("Mean Read Lengths: BA {0}bp, BCG {1}bp".format(
        BA_LENGTH, BCG_LENGTH
    ))

    BA_TOTAL = len(BA_SEQUENCES)
    BCG_TOTAL = len(BCG_SEQUENCES)
    BA_KMERS = args.ba
    BCG_KMERS = args.bcg
    LEF_KMERS = args.lef
    GENOME_SIZE = args.genome_size

    bcg_coverages = read_coverages(args.bcg_coverages)
    ba_coverages = read_coverages(args.ba_coverages)
    for bcg_coverage in bcg_coverages:
        print("Working on BCG coverage: {0}x".format(bcg_coverage))
        for ba_coverage in ba_coverages:
            print("\tWorking on BA coverage: {0}x".format(ba_coverage))
            MATCHES = 0
            start = 1
            end = args.cpu
            outputs = {}
            path = "{0}/{1}/{2}".format(args.output, bcg_coverage, ba_coverage)
            if not os.path.exists(path):
                os.makedirs(path)
            while MATCHES < args.replicates:
                with Pool(processes=args.cpu) as pool:
                    pool.map(
                        subsample,
                        [[args.working, path, bcg_coverage, ba_coverage, r]
                         for r in range(start, end + 1)]
                    )
                MATCHES = len(glob.glob("{0}/*.fasta".format(path)))
                print("\tTests: {0}, Successes: {1}.".format(end, MATCHES))
                if MATCHES < args.replicates:
                    start = end
                    end = start + args.cpu
