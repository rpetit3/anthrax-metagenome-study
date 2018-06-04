#! /usr/bin/env python3
"""Create random subsampels of input sequences."""
import argparse as ap
from multiprocessing import Pool
import os
import glob
import random
import subprocess
import numpy as np
GENOME_SIZE = None
LENGTH = None
SEQUENCES = None
TOTAL = None
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
    working_dir = opts[0]
    output_dir = opts[1]
    coverage = opts[2]
    replicate = opts[3]
    sample = opts[4]
    reads = opts[5]
    random_seed = (int(coverage * 1000000) * replicate + reads)
    basename = "{0}-{1}-{2}-{3}".format(sample, coverage, replicate,
                                        random_seed)

    if not os.path.exists('{0}/{1}-lef.txt.gz'.format(output_dir, basename)):
        fasta = []
        fasta_output = "{0}/{1}.fasta".format(working_dir, basename)
        random.seed(random_seed)
        for element in random.sample(range(TOTAL), reads):
            fasta.append(">{0}\n".format(element))
            fasta.append("{0}\n".format(SEQUENCES[element]))

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

        bcg_txt = '{0}/{1}-bcg.txt'.format(working_dir, basename)
        run_command(
            ['jellyfish', 'query', '-s', BCG_KMERS, '-o', bcg_txt, jellyfish]
        )

        lef_txt = '{0}/{1}-lef.txt'.format(working_dir, basename)
        run_command(
            ['jellyfish', 'query', '-s', LEF_KMERS, '-o', lef_txt, jellyfish]
        )

        run_command(['rm', jellyfish])
        print("\tSUCCESS: Test #: {0} Seed: {1} Reads: {2}".format(
            replicate, random_seed, reads
        ))
        run_command(['gzip', '-f', bcg_txt])
        run_command(['gzip', '-f', lef_txt])
        run_command(['gzip', '-f', ba_txt])
        run_command(['mv', '{0}.gz'.format(ba_txt), output_dir])
        run_command(['mv', '{0}.gz'.format(bcg_txt), output_dir])
        run_command(['mv', '{0}.gz'.format(lef_txt), output_dir])
    else:
        print("\tSkipping replicate: {0}, already completed".format(replicate))


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='subsample-sequences.py', conflict_handler='resolve',
        description="Create random subsamples of input sequences."
    )

    parser.add_argument('sample', type=str, metavar="SAMPLE",
                        help='Sample to subsample.')
    parser.add_argument('sequences', type=str, metavar="SEQS",
                        help='File of sequences, one per line.')
    parser.add_argument('coverages', type=str, metavar="COVERAGES",
                        help=('Coverages to subsample to.'))
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
    parser.add_argument('--length', metavar="INT", type=int, default=31,
                        help='Per line sequence length (Default 31)')
    parser.add_argument('--replicates', metavar="INT", type=int,
                        default=20,
                        help='Number of replicates per coverage (Default 100)')
    parser.add_argument('--cpu', metavar="INT", type=int, default=1,
                        help='Total number of processes to launch (Default 1)')

    args = parser.parse_args()

    SEQUENCES, LENGTH = read_sequences(args.sequences, min_length=args.length)
    print("Mean Read Length: {0}bp".format(LENGTH))

    TOTAL = len(SEQUENCES)
    BA_KMERS = args.ba
    BCG_KMERS = args.bcg
    LEF_KMERS = args.lef
    GENOME_SIZE = args.genome_size
    for coverage in read_coverages(args.coverages):
        MATCHES = 0
        start = 1
        end = args.cpu
        print("Working on coverage: {0}x ({1}, {2})".format(
            coverage, MATCHES, args.replicates
        ))
        outputs = {}
        reads = int(GENOME_SIZE * float(coverage) / LENGTH)
        path = "{0}/{1}/{2}".format(args.output, args.sample, coverage)
        if not os.path.exists(path):
            os.makedirs(path)
        with Pool(processes=args.cpu) as pool:
            pool.map(
                subsample,
                [[args.working, path, coverage, r, args.sample, reads]
                 for r in range(1, args.replicates + 1)]
            )
        MATCHES = len(glob.glob("{0}/*ba.txt.gz".format(path)))
        print("\tCompleted {0} subsamples".format(MATCHES))
