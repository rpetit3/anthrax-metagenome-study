#! /usr/bin/python3 -u
# Identify group specific kmers.
from collections import OrderedDict
import glob
import json
import logging
import os
from os.path import basename, splitext
import subprocess
import sys
import tempfile


def output_handler(output, redirect='>'):
    if output:
        return [open(output, 'w'), '{0} {1}'.format(redirect, output)]
    else:
        return [subprocess.PIPE, '']


def onfinish_handler(cmd, out, err, returncode):
    out = '\n{0}'.format(out) if out else ''
    err = '\n{0}'.format(err) if err else ''
    if returncode != 0:
        logging.error('COMMAND: {0}'.format(cmd))
        logging.error('STDOUT: {0}'.format(out))
        logging.error('STDERR: {0}'.format(err))
        logging.error('END\n'.format(err))
        raise RuntimeError(err)
    else:
        logging.info('COMMAND: {0}'.format(cmd))
        logging.info('STDOUT: {0}'.format(out))
        logging.info('STDERR: {0}'.format(err))
        logging.info('END\n'.format(err))
        return [out, err]


def byte_to_string(b):
    if b:
        return b.decode("utf-8")
    else:
        return ''


def run_command(cmd, cwd=os.getcwd(), stdout=False, stderr=False, shell=False):
    """Execute a single command and return STDOUT and STDERR."""
    stdout, stdout_str = output_handler(stdout)
    stderr, stderr_str = output_handler(stderr, redirect='2>')

    p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, cwd=cwd,
                         shell=shell)

    out, err = p.communicate()
    return onfinish_handler(
        '{0} {1} {2}'.format(' '.join(cmd), stdout_str, stderr_str),
        byte_to_string(out), byte_to_string(err), p.returncode
    )


def jellyfish_dump(sample):
    """Run Jellyfish query against a given sample."""
    cmd = ['jellyfish', 'dump', '-c', sample]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode("utf-8")
    return stdout.split('\n')


def jellyfish_query(sample, kmers):
    """Run Jellyfish query against a given sample."""
    stdout = None
    try:
        tmp = tempfile.NamedTemporaryFile()
        for kmer in kmers:
            tmp.write('>{0}\n{0}\n'.format(kmer).encode())

        p = subprocess.Popen(
            ['jellyfish', 'query', '-s', tmp.name, sample],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
        stdout = stdout.decode("utf-8")
    finally:
        tmp.close()

    return stdout.split('\n')


def parse_summary(summary):
    """Parse Summary file."""
    genomes = OrderedDict()
    cols = None
    first_anthracis = None
    with open(summary, 'r') as summary_handle:
        # Column Names:
        # accession, gi, is_bcg, is_ba, species, genome_size, description
        for line in summary_handle:
            line = line.rstrip()
            if line.startswith('#'):
                cols = line.replace('#', '').split('\t')
            else:
                row = dict(zip(cols, line.split('\t')))
                genomes[row['accession']] = row
                if not first_anthracis and row['species'] == 'anthracis':
                    first_anthracis = row['accession']
    return [genomes, first_anthracis]


def mash_order(fasta_directory, root_genome, outdir, cpu=1):
    # Sketch all of Bacillus
    if not os.path.exists('{0}/all-bacillus.msh'.format(outdir)):
        cmd = [
            'mash', 'sketch', '-k', '31', '-s', '100000', '-p', str(cpu),
            '-o', '{0}/all-bacillus'.format(outdir)
        ]
        for fasta in glob.glob('{0}/*.fasta'.format(fasta_directory)):
            cmd.append(fasta)
        run_command(cmd)

    # Mash distance between root and all other genomes
    mash_file = '{0}/root-vs-all-bacillus.dist'.format(outdir)
    if not os.path.exists(mash_file):
        run_command([
            'mash', 'dist', '-p', str(args.cpu),
            '{0}/all-bacillus.msh'.format(outdir),
            '{0}/{1}.fasta'.format(fasta_directory, root_genome)
        ], stdout=mash_file)

    distances = OrderedDict()
    cols = ['accession', 'root', 'distance', 'p-value', 'matches']
    with open(mash_file, 'r') as mash_handle:
        for line in mash_handle:
            line = line.rstrip()
            row = dict(zip(cols, line.split('\t')))
            accession = splitext(basename(row['accession']))[0]
            distances[accession] = float(row['distance'])

    sorted_distances = OrderedDict()
    for key in sorted(distances, key=lambda i: float(distances[i])):
        sorted_distances[key] = distances[key]

    return sorted_distances


def hamming_distance_filter(kmers, genomes, outdir, cpu=1):
    blastdb = "{0}/blastdb/blastdb".format(outdir)
    blast_results = "{0}/blastdb/blastn.json".format(outdir)
    run_command(['mkdir', '-p', "{0}/blastdb/".format(outdir)])
    run_command(genomes, stdout="{0}/blastdb/genomes.fasta".format(outdir))
    run_command(['makeblastdb', '-dbtype', 'nucl', '-input_type', 'fasta',
                 '-in', "{0}/blastdb/genomes.fasta".format(outdir),
                 '-out', blastdb])
    run_command([
        'blastn', '-max_hsps', '1', '-max_target_seqs', '1', '-dust', 'no',
        '-word_size', '7', '-outfmt', '15', '-query', kmers, '-db', blastdb,
        '-evalue', '10000', '-num_threads', str(cpu)
    ], stdout=blast_results)

    with open(blast_results) as blast_handle:
        json_data = json.load(blast_handle)

    hits = []
    with open("{0}/hamming-filter.txt".format(outdir), 'w') as hamming_handle:
        for entry in json_data['BlastOutput2']:
            hit = entry['report']['results']['search']
            hsp = hit['hits'][0]['hsps'][0]

            # Includes mismatches and gaps
            mismatch = hsp['align_len'] - hsp['identity']

            # Hamming distance
            hd = mismatch
            if hit['query_len'] > hsp['align_len']:
                # Include those bases that weren't aligned
                hd = hit['query_len'] - hsp['align_len'] + mismatch

            # Filter out hamming distance of 0
            if hd > 0:
                hits.append("{0}-{1}".format(hit['query_title'], hd))
            hamming_handle.write("{0}\t{1}\n".format(hit['query_title'], hd))

    return hits


def group_specific_kmers(order, organisms, fasta_directory, rrna, outdir,
                         results, is_bcg=False, cpu=1):
    info = "{0}/{1}-kmer-info.txt".format(
        results,
        'bcg' if is_bcg else 'ba'
    )
    blastdb = ['cat']
    with open(info, 'w') as info_handle:
        fasta = "{0}/{1}-kmers.fasta".format(outdir, 'bcg' if is_bcg else 'ba')
        common_kmers = None
        total = 0

        info_handle.write(
            'total\taccession\tdistance\ttotal_kmers\tgroup\tmethod\n'
        )
        print('total\taccession\tdistance\ttotal_kmers\tgroup\tmethod')
        for accession, distance in order.items():
            total += 1
            kmers = []
            group = None
            method = None
            input_file = "{0}/{1}.fasta.jf".format(fasta_directory, accession)
            if not common_kmers:
                for line in jellyfish_dump(input_file):
                    if line:
                        kmer, count = line.split(' ')
                        kmers.append(kmer)
                common_kmers = set(kmers)
            else:
                for line in jellyfish_query(input_file, list(common_kmers)):
                    if line:
                        kmer, count = line.rstrip().split(' ')

                        if int(count):
                            kmers.append(kmer)
                if organisms[accession]['is_ba'] == 'True':
                    common_kmers = common_kmers.intersection(set(kmers))
                    method = 'intersect'
                else:
                    if is_bcg and organisms[accession]['is_bcg'] == 'True':
                        common_kmers = common_kmers.intersection(set(kmers))
                        method = 'intersect'
                    else:
                        common_kmers = common_kmers.difference(set(kmers))
                        method = 'difference'
                        blastdb.append("{0}/{1}.fasta".format(
                            fasta_directory, accession
                        ))

            if organisms[accession]['is_ba'] == 'True':
                group = 'BA'
            elif organisms[accession]['is_bcg'] == 'True':
                group = 'BCG'
            else:
                group = 'Bacillus'

            k = len(common_kmers)
            info_handle.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                total, accession, distance, k, group, method
            ))
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
                total, accession, distance, k, group, method
            ))

        kmers = []
        for line in jellyfish_query(rrna, list(common_kmers)):
            if not line:
                continue
            kmer, count = line.rstrip().split(' ')
            if int(count):
                kmers.append(kmer)

        total += 1
        common_kmers = common_kmers.difference(set(kmers))
        info_handle.write("{0}\trRNA\tNA\t{1}\trRNA\tdifference\n".format(
            total, len(common_kmers)
        ))
        print("{0}\trRNA\tNA\t{1}\trRNA\tdifference".format(
            total, len(common_kmers)
        ))

        fasta = "{0}/{1}-kmers.fasta".format(outdir, 'bcg' if is_bcg else 'ba')
        with open(fasta, 'w') as fasta_handle:
            for kmer in sorted(list(common_kmers)):
                # Write as FASTA, it plays well with jellyfish query
                fasta_handle.write('>{0}\n{0}\n'.format(kmer))

        # Filter by Hamming Distance
        total += 1
        final_kmers = hamming_distance_filter(fasta, blastdb, outdir, cpu)
        info_handle.write(
            "{0}\thamming-distance\tNA\t{1}\tHD\tfilter HD =0\n".format(
                total, len(final_kmers)
            )
        )
        print("{0}\thamming-distance\tNA\t{1}\tHD\tfilter HD =0\n".format(
            total, len(final_kmers)
        ))
        fasta = "{0}/{1}-specific-kmers.fasta".format(
            results, 'bcg' if is_bcg else 'ba'
        )
        with open(fasta, 'w') as fasta_handle:
            for kmer in sorted(final_kmers):
                # Write as FASTA, it plays well with jellyfish query
                fasta_handle.write('>{0}\n{1}\n'.format(
                    kmer, kmer.split('-')[0]
                ))


if __name__ == '__main__':
    import argparse as ap
    parser = ap.ArgumentParser(
        prog='identify-group-specific-kmers.py', conflict_handler='resolve',
        description="Identify kmers only found in B. anthracis and BCG."
    )

    parser.add_argument('summary', type=str, metavar="SUMMARY",
                        help='Summary of Bacillus genomes.')
    parser.add_argument('fasta', type=str, metavar="FASTA",
                        help='Directory of completed genomes in FASTA format')
    parser.add_argument('rrna', type=str, metavar="RRNA",
                        help=('A dump of rRNA kmers.'))
    parser.add_argument('directory', type=str, metavar="OUTPUT",
                        help='Directory for analysis results.')
    parser.add_argument('results', type=str, metavar="OUTPUT",
                        help='Directory to output final results.')
    parser.add_argument('--cpu', metavar="INT", type=int, default=1,
                        help='Number of processors to use (Default 1)')

    args = parser.parse_args()
    # Setup logs
    logging.basicConfig(filename='identify-group-specific-kmers.txt',
                        filemode='w', level=logging.INFO)
    genomes, first_anthracis = parse_summary(args.summary)

    # Run mash to get genome order, starting with first B. anthracis genome in
    # the summary file
    genome_order = mash_order(args.fasta, first_anthracis, args.directory,
                              cpu=args.cpu)

    # B. anthracis specific kmers
    group_specific_kmers(genome_order, genomes, args.fasta, args.rrna,
                         args.directory, args.results, cpu=args.cpu)
    # BCG specific kmers
    group_specific_kmers(genome_order, genomes, args.fasta, args.rrna,
                         args.directory, args.results, is_bcg=True,
                         cpu=args.cpu)
