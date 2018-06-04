#! /usr/bin/env python3
"""Parse through the simulated sequencing group specific kmer counts."""
import argparse as ap
from collections import OrderedDict
import glob
import gzip
import os
import sys
import time
import numpy as np
import multiprocessing as mp

SAMPLES = OrderedDict()
KMERS = {}
HAMMING = OrderedDict()

SAMPLE_COLS = [
    'sample', 'is_bcg', 'is_ba', 'has_lethal', 'simulated_coverage', 'group',
    'total_kmers', 'tp', 'tn', 'fp', 'fn',
    'kmer_cov_min', 'kmer_cov_mean', 'kmer_cov_median', 'kmer_cov_max',
    'non_zero_kmer_cov_min', 'non_zero_kmer_cov_mean',
    'non_zero_kmer_cov_median', 'non_zero_kmer_cov_max'
]

KMER_COLS = [
    'kmer', 'simulated_coverage', 'group', 'hamming_distance',
    'tp', 'tn', 'fp', 'fn',
    'group_kmer_cov_min',
    'group_kmer_cov_mean',
    'group_kmer_cov_median',
    'group_kmer_cov_max',
    'non_zero_group_kmer_cov_min',
    'non_zero_group_kmer_cov_mean',
    'non_zero_group_kmer_cov_median',
    'non_zero_group_kmer_cov_max',
    'outgroup_kmer_cov_min',
    'outgroup_kmer_cov_mean',
    'outgroup_kmer_cov_median',
    'outgroup_kmer_cov_max',
    'non_zero_outgroup_kmer_cov_min',
    'non_zero_outgroup_kmer_cov_mean',
    'non_zero_outgroup_kmer_cov_median',
    'non_zero_outgroup_kmer_cov_max'
]


def get_group_status(sample, group):
    """Return if a sample is within a group or not."""
    within_group = None
    if group == 'ba':
        within_group = True if SAMPLES[sample]['is_ba'] == 'True' else False
    elif group == 'bcg':
        within_group = True if SAMPLES[sample]['is_bcg'] == 'True' else False
    else:
        # lef
        within_group = True if SAMPLES[sample]['has_lethal'] else False

    return within_group


def get_coverage_stats(coverage):
    """Return summary stats of a set of coverages."""
    non_zero = [c for c in coverage if c]
    np_array = np.array(coverage)
    non_zero_array = np.array(non_zero)
    return {
        'min': min(coverage) if coverage else 0,
        'median': int(np.median(np_array)) if coverage else 0,
        'mean': "{0:.4f}".format(np.mean(np_array)) if coverage else 0,
        'max': max(coverage) if coverage else 0,
        'non_zero_min': min(non_zero_array) if non_zero else 0,
        'non_zero_median': int(np.median(non_zero_array)) if non_zero else 0,
        'non_zero_mean': int(round(np.mean(non_zero_array))) if non_zero else 0,
        'non_zero_max': max(non_zero_array) if non_zero else 0,
    }


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    return ''.join([complement[b] for b in seq[::-1]])


def parse_counts(counts, sample, coverage, group, skip_kmers=False,
                 filter_kmers=False):
    """Parse kmer counts."""
    within_group = get_group_status(sample, group)
    sample_row = {'coverages': [], 'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
    with gzip.open(counts, 'r') as count_handle:
        for line in count_handle:
            kmer, count = line.decode().rstrip().split()
            count = int(count)
            parse = True

            if filter_kmers:
                parse = kmer in KMERS or reverse_complement(kmer) in KMERS
            elif not skip_kmers:
                if kmer not in KMERS:
                    kmer = reverse_complement(kmer)
                if within_group:
                    KMERS[kmer][coverage]['group_coverages'].append(count)
                    if count:
                        KMERS[kmer][coverage]['tp'] += 1
                    else:
                        KMERS[kmer][coverage]['fn'] += 1
                else:
                    KMERS[kmer][coverage]['outgroup_coverages'].append(count)
                    if count:
                        KMERS[kmer][coverage]['fp'] += 1
                    else:
                        KMERS[kmer][coverage]['tn'] += 1

            if parse:
                sample_row['coverages'].append(count)
                if within_group:
                    if count:
                        sample_row['tp'] += 1
                    else:
                        sample_row['fn'] += 1
                else:
                    if count:
                        sample_row['fp'] += 1
                    else:
                        sample_row['tn'] += 1

    coverage_stats = get_coverage_stats(sample_row['coverages'])
    SAMPLES[sample]['results'].append({
        'simulated_coverage': coverage,
        'within_group': within_group,
        'tp': sample_row['tp'],
        'tn': sample_row['tn'],
        'fp': sample_row['fp'],
        'fn': sample_row['fn'],
        'kmer_cov_min': coverage_stats['min'],
        'kmer_cov_mean': coverage_stats['mean'],
        'kmer_cov_median': coverage_stats['median'],
        'kmer_cov_max': coverage_stats['max'],
        'non_zero_kmer_cov_min': coverage_stats['non_zero_min'],
        'non_zero_kmer_cov_mean': coverage_stats['non_zero_mean'],
        'non_zero_kmer_cov_median': coverage_stats['non_zero_median'],
        'non_zero_kmer_cov_max': coverage_stats['non_zero_max'],
    })


def parse_kmers(kmers, coverages, skip_kmers=False, has_hamming=True):
    with open(kmers, 'r') as kmer_handle:
        for line in kmer_handle:
            if line.startswith(">"):
                line = line.rstrip().replace(">", "")
                kmer, distance = line.split("-")
                if not has_hamming:
                    distance = False
                KMERS[kmer] = OrderedDict()
                HAMMING[kmer] = distance
                if not skip_kmers:
                    for coverage in coverages:
                        KMERS[kmer][coverage] = {
                            'group_coverages': [], 'outgroup_coverages': [],
                            'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0
                        }


def parse_summary(summary):
    """Parse Summary file."""
    cols = None
    with open(summary, 'r') as summary_handle:
        # Column Names:
        # accession, gi, is_bcg, is_ba, species, genome_size, description
        for line in summary_handle:
            line = line.rstrip()
            if line.startswith('#'):
                cols = line.replace('#', '').split('\t')
            else:
                row = dict(zip(cols, line.split('\t')))
                SAMPLES[row['accession']] = row
                if row['accession'] == 'NZ_CP009941':
                    # NZ_CP009941 - Bacillus cereus  w/ lef on chromosome
                    SAMPLES[row['accession']]['has_lethal'] = True
                else:
                    SAMPLES[row['accession']]['has_lethal'] = False
                SAMPLES[row['accession']]['results'] = []


def print_sample_summary(file_output):
    """Print the final per sample summaries."""
    with open(file_output, 'w') as output_handle:
        output_handle.write(("\t".join(SAMPLE_COLS)))
        output_handle.write("\n")
        for sample in SAMPLES:
            if SAMPLES[sample]['results']:
                for result in SAMPLES[sample]['results']:
                    row = {
                        'sample': sample,
                        'is_bcg': SAMPLES[sample]['is_bcg'],
                        'is_ba': SAMPLES[sample]['is_ba'],
                        'has_lethal': SAMPLES[sample]['has_lethal'],
                        'simulated_coverage': result['simulated_coverage'],
                        'group': args.group,
                        'within_group': result['within_group'],
                        'total_kmers': total_kmers,
                        'tp': result['tp'],
                        'tn': result['tn'],
                        'fp': result['fp'],
                        'fn': result['fn'],
                        'kmer_cov_min': result['kmer_cov_min'],
                        'kmer_cov_mean': result['kmer_cov_mean'],
                        'kmer_cov_median': result['kmer_cov_median'],
                        'kmer_cov_max': result['kmer_cov_max'],
                        'non_zero_kmer_cov_min': result['non_zero_kmer_cov_min'],
                        'non_zero_kmer_cov_mean': result['non_zero_kmer_cov_mean'],
                        'non_zero_kmer_cov_median': result['non_zero_kmer_cov_median'],
                        'non_zero_kmer_cov_max': result['non_zero_kmer_cov_max']
                    }
                    output_handle.write(("\t".join([
                        str(row[col]) for col in SAMPLE_COLS
                    ])))
                    output_handle.write("\n")


def print_kmer_summary(file_output):
    """Print the final per kmer summaries."""
    with open(file_output, 'w') as output_handle:
        output_handle.write(("\t".join(KMER_COLS)))
        output_handle.write("\n")
        for kmer, coverages in KMERS.items():
            for coverage in coverages:
                within_group = get_coverage_stats(
                    KMERS[kmer][coverage]['group_coverages']
                )
                outgroup = get_coverage_stats(
                    KMERS[kmer][coverage]['outgroup_coverages']
                )
                row = {
                    'kmer': kmer,
                    'simulated_coverage': coverage,
                    'group': args.group,
                    'hamming_distance': HAMMING[kmer],
                    'tp': KMERS[kmer][coverage]['tp'],
                    'tn': KMERS[kmer][coverage]['tn'],
                    'fp': KMERS[kmer][coverage]['fp'],
                    'fn': KMERS[kmer][coverage]['fn'],
                    'group_kmer_cov_min': within_group['min'],
                    'group_kmer_cov_mean': within_group['mean'],
                    'group_kmer_cov_median': within_group['median'],
                    'group_kmer_cov_max': within_group['max'],

                    'non_zero_group_kmer_cov_min': within_group['non_zero_min'],
                    'non_zero_group_kmer_cov_mean': within_group['non_zero_mean'],
                    'non_zero_group_kmer_cov_median': within_group['non_zero_median'],
                    'non_zero_group_kmer_cov_max': within_group['non_zero_max'],

                    'outgroup_kmer_cov_min': outgroup['min'],
                    'outgroup_kmer_cov_mean': outgroup['mean'],
                    'outgroup_kmer_cov_median': outgroup['median'],
                    'outgroup_kmer_cov_max': outgroup['max'],

                    'non_zero_outgroup_kmer_cov_min': outgroup['non_zero_min'],
                    'non_zero_outgroup_kmer_cov_mean': outgroup['non_zero_mean'],
                    'non_zero_outgroup_kmer_cov_median': outgroup['non_zero_median'],
                    'non_zero_outgroup_kmer_cov_max': outgroup['non_zero_max'],
                }
                output_handle.write(("\t".join([
                    str(row[col]) for col in KMER_COLS
                ])))
                output_handle.write("\n")


def read_lines(input_file):
    """Return lines in a text file as a list."""
    lines = []
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            lines.append(line.rstrip())
    return lines


def parse_filter_kmers(kmers):
    with open(kmers, 'r') as kmer_handle:
        for line in kmer_handle:
            if line.startswith(">"):
                line = line.rstrip().replace(">", "")
                KMERS[line.split("-")[0]] = True


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='summarize-kmer-counts.py', conflict_handler='resolve',
        description=("Summarize kmer counts of each simulation.")
    )

    parser.add_argument('summary', type=str, metavar="SUMMARY",
                        help='Summary of Bacillus genomes.')
    parser.add_argument('directory', type=str, metavar="SIMUALTION_DIR",
                        help='Directory with group specific 31-mer counts.')
    parser.add_argument('group', type=str, metavar="GROUP",
                        help='Which group to parse (ba, bcg or lef).')
    parser.add_argument('kmers', type=str, metavar="KMERS",
                        help='Group specific k-mers.')
    parser.add_argument('coverages', type=str, metavar="COVERAGES",
                        help=('Coverages to subsample to.'))
    parser.add_argument('outdir', type=str, metavar="OUTDIR",
                        help='Directory to output to.')
    parser.add_argument('--cpu', default=1, type=int, metavar="INT",
                        help='Number of cores to use (Default: 1)')
    parser.add_argument('--single_sample', type=str, metavar="STR",
                        help='Process a single sample.')
    parser.add_argument('--skip_kmers', action='store_true', default=False,
                        help='Skip kmer processing.')
    parser.add_argument('--filter', action='store_true', default=False,
                        help='Filter counts based on input kmers.')
    args = parser.parse_args()

    if args.group not in ['ba', 'bcg', 'lef']:
        raise Exception("GROUPS must be 'ba', 'bcg' or 'lef'")

    coverages = read_lines(args.coverages)
    print("Parsing Summary")
    parse_summary(args.summary)
    print("Parsing Kmers")
    if args.filter:
        print("Filtering Kmers")
        args.skip_kmers = True
        parse_filter_kmers(args.kmers)
    else:
        print("Parsing Kmers")
        parse_kmers(args.kmers, coverages, skip_kmers=args.skip_kmers,
                    has_hamming=False if args.group == 'lef' else True)
    total_kmers = len(KMERS)

    current = 1
    samples = list(SAMPLES.keys())
    if args.single_sample:
        samples = [args.single_sample]
    total = len(samples)
    for sample in samples:
        path = "{0}/{1}".format(args.directory, sample)
        if os.path.exists(path):
            print("Working on {0} ({1} of {2})".format(sample, current, total))
            current += 1
            count_files = sorted(glob.glob(
                "{0}/*-{1}.txt.gz".format(path, args.group)
            ))
            for count_file in count_files:
                coverage = os.path.basename(count_file).split('-')[1]
                parse_counts(count_file, sample, coverage, args.group,
                             skip_kmers=args.skip_kmers,
                             filter_kmers=args.filter)

    print("Output sample summary")
    if args.single_sample:
        print_sample_summary("{0}/count-summary-{1}-{2}.txt".format(
            args.outdir, args.single_sample, args.group
        ))
    else:
        print_sample_summary("{0}/count-summary-sample-{1}.txt".format(
            args.outdir, args.group
        ))

    if not args.skip_kmers:
        print("Output kmer summary")
        if args.single_sample:
            print_kmer_summary("{0}/count-summary-kmer-{1}-{2}.txt".format(
                args.outdir, args.single_sample, args.group
            ))
        else:
            print_kmer_summary("{0}/count-summary-kmer-{1}.txt".format(
                args.outdir, args.group
            ))
