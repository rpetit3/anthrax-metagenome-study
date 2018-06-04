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
KMERS = OrderedDict()
HAMMING = OrderedDict()

SAMPLE_COLS = [
    'run', 'is_bcg', 'is_ba', 'has_lethal', 'group',
    'total_kmers', 'tp', 'tn', 'fp', 'fn',
    'kmer_cov_min', 'kmer_cov_mean', 'kmer_cov_median', 'kmer_cov_max',
    'non_zero_kmer_cov_min', 'non_zero_kmer_cov_mean',
    'non_zero_kmer_cov_median', 'non_zero_kmer_cov_max'
]

META_SAMPLE_COLS = [
    'run', 'is_bcg', 'is_ba', 'has_lethal', 'group', 'total_kmers',
    'hit', 'miss', 'kmer_cov_min', 'kmer_cov_mean', 'kmer_cov_median',
    'kmer_cov_max', 'non_zero_kmer_cov_min', 'non_zero_kmer_cov_mean',
    'non_zero_kmer_cov_median', 'non_zero_kmer_cov_max'
]

KMER_COLS = [
    'kmer', 'group', 'hamming_distance', 'tp', 'tn', 'fp', 'fn',
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

META_KMER_COLS = [
    'kmer', 'group', 'hamming_distance', 'total', 'hit', 'miss',
    'kmer_cov_min', 'kmer_cov_mean', 'kmer_cov_median', 'kmer_cov_max',
    'non_zero_kmer_cov_min', 'non_zero_kmer_cov_mean',
    'non_zero_kmer_cov_median', 'non_zero_kmer_cov_max'
]


def get_group_status(sample, group):
    """Return if a sample is within a group or not."""
    within_group = None
    if group == 'ba':
        within_group = True if SAMPLES[sample]['is_ba'] else False
    elif group == 'bcg':
        within_group = True if SAMPLES[sample]['is_bcg'] else False
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
        'mean': np.mean(np_array) if coverage else 0,
        'max': max(coverage) if coverage else 0,
        'non_zero_min': min(non_zero_array) if non_zero else 0,
        'non_zero_median': int(np.median(non_zero_array)) if non_zero else 0,
        'non_zero_mean': np.mean(non_zero_array) if non_zero else 0,
        'non_zero_max': max(non_zero_array) if non_zero else 0,
    }


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    return ''.join([complement[b] for b in seq[::-1]])


def parse_counts(counts, sample, group, is_metagenome=False,
                 exclude_singletons=False):
    """Parse kmer counts."""
    within_group = get_group_status(sample, group)
    sample_row = {'coverages': [], 'total_kmers': 0, 'tp': 0, 'tn': 0,
                  'fp': 0, 'fn': 0}
    if is_metagenome:
        sample_row = {'coverages': [], 'total_kmers': 0, 'hit': 0, 'miss': 0}
    with gzip.open(counts, 'r') as count_handle:
        for line in count_handle:
            kmer, count = line.decode().rstrip().split()

            if kmer not in KMERS:
                kmer = reverse_complement(kmer)

            if kmer in KMERS:
                count = int(count)
                if exclude_singletons and count == 1:
                    continue

                sample_row['coverages'].append(count)
                sample_row['total_kmers'] += 1
                if is_metagenome:
                    KMERS[kmer]['coverages'].append(count)
                    KMERS[kmer]['total'] += count
                    if count:
                        sample_row['hit'] += 1
                        KMERS[kmer]['hit'] += 1
                    else:
                        sample_row['miss'] += 1
                        KMERS[kmer]['miss'] += 1
                else:
                    if within_group:
                        KMERS[kmer]['group_coverages'].append(count)
                        if count:
                            sample_row['tp'] += 1
                            KMERS[kmer]['tp'] += 1
                        else:
                            sample_row['fn'] += 1
                            KMERS[kmer]['fn'] += 1
                    else:
                        KMERS[kmer]['outgroup_coverages'].append(count)
                        if count:
                            sample_row['fp'] += 1
                            KMERS[kmer]['fp'] += 1
                        else:
                            sample_row['tn'] += 1
                            KMERS[kmer]['tn'] += 1

    coverage_stats = get_coverage_stats(sample_row['coverages'])
    result = {
        'within_group': within_group,
        'total_kmers': sample_row['total_kmers'],
        'kmer_cov_min': coverage_stats['min'],
        'kmer_cov_mean': coverage_stats['mean'],
        'kmer_cov_median': coverage_stats['median'],
        'kmer_cov_max': coverage_stats['max'],
        'non_zero_kmer_cov_min': coverage_stats['non_zero_min'],
        'non_zero_kmer_cov_mean': coverage_stats['non_zero_mean'],
        'non_zero_kmer_cov_median': coverage_stats['non_zero_median'],
        'non_zero_kmer_cov_max': coverage_stats['non_zero_max'],
    }
    if is_metagenome:
        result['hit'] = sample_row['hit']
        result['miss'] = sample_row['miss']
    else:
        result['tp'] = sample_row['tp']
        result['tn'] = sample_row['tn']
        result['fp'] = sample_row['fp']
        result['fn'] = sample_row['fn']

    SAMPLES[sample]['results'].append(result)


def parse_kmers(kmers, has_hamming=True, is_metagenome=False):
    with open(kmers, 'r') as kmer_handle:
        for line in kmer_handle:
            if line.startswith(">"):
                line = line.rstrip().replace(">", "")
                kmer, distance = line.split("-")
                if not has_hamming:
                    distance = False
                KMERS[kmer] = OrderedDict()
                if is_metagenome:
                    KMERS[kmer] = {
                        'distance': distance,
                        'coverages': [],
                        'total': 0,
                        'hit': 0,
                        'miss': 0,
                    }
                else:
                    KMERS[kmer] = {
                        'distance': distance,
                        'group_coverages': [],
                        'outgroup_coverages': [],
                        'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0
                    }


def parse_summary(summary, is_bcg=False, is_metagenome=False):
    """Parse Summary file."""
    with open(summary, 'r') as summary_handle:
        # Column Names:
        # run_accession, organism
        for line in summary_handle:
            if not line.startswith('#'):
                line = line.rstrip()
                if is_metagenome:
                    SAMPLES[line] = {
                        'is_bcg': False,
                        'is_ba': False,
                        'has_lethal': False,
                        'results': []
                    }
                else:
                    accession, organism = line.split('\t')
                    SAMPLES[accession] = {
                        'is_bcg': is_bcg,
                        'is_ba': False,
                        'has_lethal': False,
                        'results': []
                    }
                    if organism.startswith("Bacillus anthracis"):
                        SAMPLES[accession]['is_ba'] = True
                        SAMPLES[accession]['has_lethal'] = True
                    elif "biovar anthracis" in organism:
                        SAMPLES[accession]['has_lethal'] = True


def print_sample_summary(file_output, group, is_metagenome=False):
    """Print the final per sample summaries."""
    with open(file_output, 'w') as output_handle:
        if is_metagenome:
            output_handle.write(("\t".join(META_SAMPLE_COLS)))
            output_handle.write("\n")
        else:
            output_handle.write(("\t".join(SAMPLE_COLS)))
            output_handle.write("\n")
        for sample in SAMPLES:
            if SAMPLES[sample]['results']:
                for result in SAMPLES[sample]['results']:
                    row = None
                    if is_metagenome:
                        row = {
                            'run': sample,
                            'is_bcg': SAMPLES[sample]['is_bcg'],
                            'is_ba': SAMPLES[sample]['is_ba'],
                            'has_lethal': SAMPLES[sample]['has_lethal'],
                            'group': group,
                            'within_group': result['within_group'],
                            'total_kmers': result['total_kmers'],
                            'hit': result['hit'],
                            'miss': result['miss'],
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
                            str(row[col]) for col in META_SAMPLE_COLS
                        ])))
                        output_handle.write("\n")
                    else:
                        row = {
                            'run': sample,
                            'is_bcg': SAMPLES[sample]['is_bcg'],
                            'is_ba': SAMPLES[sample]['is_ba'],
                            'has_lethal': SAMPLES[sample]['has_lethal'],
                            'group': group,
                            'within_group': result['within_group'],
                            'total_kmers': result['total_kmers'],
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


def print_kmer_summary(file_output, group, is_metagenome=False):
    """Print the final per kmer summaries."""
    with open(file_output, 'w') as output_handle:
        if is_metagenome:
            output_handle.write(("\t".join(META_KMER_COLS)))
            output_handle.write("\n")
        else:
            output_handle.write(("\t".join(KMER_COLS)))
            output_handle.write("\n")
        for kmer, vals in KMERS.items():
            if is_metagenome:
                coverages = get_coverage_stats(vals['coverages'])
                row = {
                    'kmer': kmer,
                    'group': group,
                    'hamming_distance': vals['distance'],
                    'total': vals['total'],
                    'hit': vals['hit'],
                    'miss': vals['miss'],
                    'kmer_cov_min': coverages['min'],
                    'kmer_cov_mean': coverages['mean'],
                    'kmer_cov_median': coverages['median'],
                    'kmer_cov_max': coverages['max'],
                    'non_zero_kmer_cov_min': coverages['non_zero_min'],
                    'non_zero_kmer_cov_mean': coverages['non_zero_mean'],
                    'non_zero_kmer_cov_median': coverages['non_zero_median'],
                    'non_zero_kmer_cov_max': coverages['non_zero_max']
                }
                output_handle.write(("\t".join([
                    str(row[col]) for col in META_KMER_COLS
                ])))
                output_handle.write("\n")
            else:
                within_group = get_coverage_stats(vals['group_coverages'])
                outgroup = get_coverage_stats(vals['outgroup_coverages'])
                row = {
                    'kmer': kmer,
                    'group': group,
                    'hamming_distance': vals['distance'],
                    'tp': vals['tp'],
                    'tn': vals['tn'],
                    'fp': vals['fp'],
                    'fn': vals['fn'],
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


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='summarize-public-data-counts.py', conflict_handler='resolve',
        description=("Summarize kmer counts of of public data sequencing.")
    )

    parser.add_argument('summary', type=str, metavar="SUMMARY",
                        help='Summary of runs accessions.')
    parser.add_argument('directory', type=str, metavar="SIMUALTION_DIR",
                        help='Directory with group specific 31-mer counts.')
    parser.add_argument('group', type=str, metavar="GROUP",
                        help='Which group to parse (ba, bcg or lef).')
    parser.add_argument('kmers', type=str, metavar="KMERS",
                        help='Group specific k-mers.')
    parser.add_argument('prefix', type=str, metavar="PREFIX",
                        help='Prefix for saving output.')
    parser.add_argument('--is_bcg', action='store_true', default=False,
                        help='Samples are from Bacillus cereus group.')
    parser.add_argument('--is_metagenome', action='store_true', default=False,
                        help='Samples are metegenomic sequencing.')
    parser.add_argument('--exclude_singletons', action='store_true',
                        default=False, help='Skip kmers counted only once.')
    args = parser.parse_args()

    if args.group not in ['ba', 'bcg', 'lef']:
        raise Exception("GROUPS must be 'ba', 'bcg' or 'lef'")

    print("Parsing Summary")
    parse_summary(args.summary, is_bcg=args.is_bcg,
                  is_metagenome=args.is_metagenome)

    print("Parsing Kmers")
    parse_kmers(args.kmers, has_hamming=False if args.group == 'lef' else True,
                is_metagenome=args.is_metagenome)

    current = 1
    total = len(SAMPLES)
    for sample in SAMPLES:
        print("Working on {0} ({1} of {2})".format(sample, current, total))
        current += 1
        count_file = "{0}/{1}-{2}.txt.gz".format(
            args.directory, sample, args.group
        )

        if os.path.exists(count_file):
            parse_counts(count_file, sample, args.group,
                         is_metagenome=args.is_metagenome,
                         exclude_singletons=args.exclude_singletons)

    print("Output sample summary")
    print_sample_summary("{0}-count-summary-sample-{1}.txt".format(
        args.prefix, args.group
    ), args.group, is_metagenome=args.is_metagenome)
    print("Output kmer summary")
    print_kmer_summary("{0}-count-summary-kmer-{1}.txt".format(
        args.prefix, args.group
    ), args.group, is_metagenome=args.is_metagenome)
