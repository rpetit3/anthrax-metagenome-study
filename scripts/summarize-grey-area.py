#! /usr/bin/env python3
"""Parse through the subsample counts of group specific kmer counts."""
import argparse as ap
from collections import OrderedDict
import glob
import gzip
import numpy as np
import os

SAMPLES = {
    'ba': OrderedDict(),
    'bcg': OrderedDict(),
    'lef': OrderedDict()
}

SAMPLE_COLS = [
    'replicate', 'group', 'bcg_coverage', 'ba_coverage', 'is_bcg', 'is_ba',
    'has_lethal', 'total_kmers', 'tp', 'tn', 'fp', 'fn',
    'kmer_cov_min', 'kmer_cov_mean', 'kmer_cov_median', 'kmer_cov_max',
    'non_zero_kmer_cov_min', 'non_zero_kmer_cov_mean',
    'non_zero_kmer_cov_median', 'non_zero_kmer_cov_max'
]

KMERS = {}


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


def parse_counts(counts, sample, group, bcg_coverage, ba_coverage):
    """Parse kmer counts."""
    sample_row = {
        'coverages': [], 'total_kmers': 0, 'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0
    }

    cols = []
    with gzip.open(counts, 'r') as count_handle:
        for line in count_handle:
            kmer, count = line.decode().rstrip().split()
            count = int(count)
            parse = False
            if KMERS and group == 'ba':
                if kmer in KMERS or reverse_complement(kmer) in KMERS:
                    parse = True
            else:
                parse = True

            if parse:
                sample_row['coverages'].append(count)
                sample_row['total_kmers'] += 1
                if ba_coverage:
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
    SAMPLES[group][sample] = {
        'replicate': sample,
        'group': group,
        'bcg_coverage': bcg_coverage,
        'ba_coverage': ba_coverage,
        'is_bcg': True,
        'is_ba': True if ba_coverage else False,
        'has_lethal': True if ba_coverage else False,
        'total_kmers': sample_row['total_kmers'],
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
    }


def print_sample_summary(file_output, group):
    """Print the final per sample summaries."""
    with open(file_output, 'w') as output_handle:
        output_handle.write(("\t".join(SAMPLE_COLS)))
        output_handle.write("\n")
        for sample, vals in SAMPLES[group].items():
            output_handle.write(("\t".join([
                str(vals[col]) for col in SAMPLE_COLS
            ])))
            output_handle.write("\n")


def parse_kmers(kmers):
    with open(kmers, 'r') as kmer_handle:
        for line in kmer_handle:
            if line.startswith(">"):
                line = line.rstrip().replace(">", "")
                KMERS[line.split("-")[0]] = True


def parse_directory(path):
    """"Pull out BCG and BA subsample coverages."""
    folders = path.replace("/", " ").rstrip().split()
    return [float(folders[-2]), float(folders[-1])]


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='summarize-grey-area.py', conflict_handler='resolve',
        description=("Summarize kmer counts of of grey area subsampling.")
    )

    parser.add_argument('counts', type=str, metavar="COUNTS",
                        help='Directory of counts.')
    parser.add_argument('output', type=str, metavar="OUTPUT_DIR",
                        help='Directory to output summaries.')
    parser.add_argument('--filter', type=str, metavar="OUTPUT_DIR",
                        help='Directory to output summaries.')
    args = parser.parse_args()

    if args.filter:
        parse_kmers(args.filter)

    bcg_coverage, ba_coverage = parse_directory(args.counts)
    if ba_coverage <= 5.0:
        for group in ['ba', 'bcg', 'lef']:
            glob_string = "{0}/*-{1}.txt.gz".format(args.counts, group)
            for counts in glob.glob(glob_string):
                replicate = os.path.basename(counts).split("-")[1]
                parse_counts(counts, replicate, group, bcg_coverage,
                             ba_coverage)

            output_path = "{0}/{1}".format(args.output, group)
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            print_sample_summary("{0}/{1}-{2}-{3}.txt".format(
                output_path, bcg_coverage, ba_coverage, group
            ), group)
