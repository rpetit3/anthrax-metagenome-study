#! /usr/bin/env python3
"""Create random subsampels of input sequences."""
import argparse as ap
from multiprocessing import Pool
import os
import glob
import subprocess
GENOME_SIZE = None
LENGTH = None
SEQUENCES = None
TOTAL = None
BA_KMERS = None
BCG_KMERS = None
LEF_KMERS = None
MATCHES = 0


def read_coverages(input_file):
    """Return coverage values."""
    coverages = []
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            coverages.append(float(line.rstrip()))
    return coverages


def has_hit(counts):
    with open(counts, 'r') as count_handle:
        for line in count_handle:
            if int(line.rstrip().split()[1]):
                return True
    return False


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
    input_fasta = opts[4]
    sample = opts[5]
    random_seed = (int(coverage * 10000) * replicate)
    basename = "{0}-{1}-{2}-{3}".format(sample, coverage, replicate,
                                        random_seed)
    if not os.path.exists('{0}/{1}-lef.txt.gz'.format(output_dir, basename)):
        # Count kmers
        art_prefix = "{0}/{1}".format(working_dir, basename)
        art_fastq = "{0}/{1}.fq".format(working_dir, basename)
        run_command(['art_illumina', '-l', '100', '-f', str(coverage), '-na',
                     '-ss', 'HS20', '-rs', str(random_seed), '-i', input_fasta,
                     '-o', art_prefix])

        jellyfish = '{0}/{1}.jf'.format(working_dir, basename)
        run_command(['jellyfish', 'count', '-C', '-t', '4', '-m', '31',
                     '-s', '5M', '-o', jellyfish, art_fastq])
        run_command(['rm', art_fastq])
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
            print("\tSUCCESS: Test #: {0} Seed: {1}".format(
                replicate, random_seed
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
    else:
        print("\tSkipping replicate: {0}, already completed".format(replicate))


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='simpulate-single-sample.py', conflict_handler='resolve',
        description="Create random subsamples of input sequences."
    )

    parser.add_argument('sample', type=str, metavar="SAMPLE",
                        help='Sample to subsample.')
    parser.add_argument('fasta', type=str, metavar="FASTA",
                        help='FASTA sequence to subsample.')
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
    parser.add_argument('--replicates', metavar="INT", type=int,
                        default=100,
                        help='Number of replicates per coverage (Default 100)')
    parser.add_argument('--cpu', metavar="INT", type=int, default=1,
                        help='Total number of processes to launch (Default 1)')

    args = parser.parse_args()
    BA_KMERS = args.ba
    BCG_KMERS = args.bcg
    LEF_KMERS = args.lef

    coverages = read_coverages(args.coverages)
    tests = []
    for coverage in coverages:
        print("Working on coverage: {0}x".format(coverage))
        path = "{0}/{1}/{2}".format(args.output, args.sample, coverage)
        if not os.path.exists(path):
            os.makedirs(path)
        MATCHES = 0
        start = 1
        end = args.replicates
        while MATCHES < args.replicates:
            with Pool(processes=args.cpu) as pool:
                pool.map(
                    subsample,
                    [[args.working, path, coverage, r, args.fasta, args.sample]
                     for r in range(start, end + 1)]
                )
            MATCHES = len(glob.glob("{0}/*ba.txt.gz".format(path)))

            print("\tTests: {0}, Successes: {1}.".format(end, MATCHES))
            if MATCHES < args.replicates:
                start = end
                end = start + args.replicates
        tests.append([str(coverage), str(end), str(MATCHES)])
    with open('{0}-summary.txt'.format(args.sample), 'w') as output_handle:
        output_handle.write("coverage\ttests\tsuccesses\n")
        for row in tests:
            output_handle.write("{0}\n".format("\t".join(row)))
