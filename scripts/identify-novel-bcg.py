#! /usr/bin/env python3
# Identify nonBCG members that should be apart of BCG.
from collections import OrderedDict
import glob
import logging
import os
from os.path import basename, splitext
import subprocess
import sys


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


if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='identify-novel-bcg.py',
        conflict_handler='resolve',
        description=('Identify nonBCG members that should be apart of BCG.'))

    parser.add_argument('summary', type=str, metavar="SUMMARY",
                        help='Summary of Bacillus genomes.')
    parser.add_argument('directory', type=str, metavar="OUTPUT",
                        help='Directory to output results.')
    parser.add_argument('fasta', type=str, metavar="FASTA",
                        help='Directory of completed genomes in FASTA format')
    parser.add_argument('--cpu', metavar="INT", type=int, default=1,
                        help='Number of processors to use (Default 1)')
    args = parser.parse_args()
    # Setup logs
    logging.basicConfig(filename='identify-novel-bcg.txt',
                        filemode='w', level=logging.INFO)

    # Make tempory directories for Bacillus and BCG
    fasta_dir = '{0}/{1}'.format(os.getcwd(), args.fasta)
    outdir = '{0}/{1}'.format(os.getcwd(), args.directory)
    bcg_fasta = "{0}/bcg".format(outdir)
    bacillus_fasta = "{0}/bacillus".format(outdir)
    run_command(['mkdir', '-p', bcg_fasta])
    run_command(['mkdir', '-p', bacillus_fasta])

    # Parse Summary files and create symbolic links
    genomes = OrderedDict()
    cols = None
    with open(args.summary, 'r') as summary_handle:
        # Column Names:
        # accession, gi, is_bcg, is_ba, species, genome_size, description
        for line in summary_handle:
            line = line.rstrip()
            if line.startswith('#'):
                cols = line.replace('#', '').split('\t')
            else:
                row = dict(zip(cols, line.split('\t')))
                genomes[row['accession']] = row
                source = "{0}/{1}.fasta".format(fasta_dir, row['accession'])
                target = "{0}/{1}.fasta".format(
                    bcg_fasta if row['is_bcg'] == 'True' else bacillus_fasta,
                    row['accession']
                )
                run_command(['cp', source, target])

    # Mash BCG against BCG
    if not os.path.exists('{0}/bcg.msh'.format(outdir)):
        cmd = [
            'mash', 'sketch', '-k', '31', '-s', '10000', '-p', str(args.cpu),
            '-o', '{0}/bcg'.format(outdir)
        ]
        for fasta in glob.glob('{}/*.fasta'.format(bcg_fasta)):
            cmd.append(fasta)
        run_command(cmd)

    if not os.path.exists('{0}/bcg.dist'.format(outdir)):
        run_command([
            'mash', 'dist', '-p', str(args.cpu), '{0}/bcg.msh'.format(outdir),
            '{0}/bcg.msh'.format(outdir)
        ], stdout='{0}/bcg.dist'.format(outdir))

    # Read the within BCG member distances
    within_bcg_distances = []
    with open('{0}/bcg.dist'.format(outdir), 'r') as distance_handle:
        for line in distance_handle:
            within_bcg_distances.append(float(line.split('\t')[2]))

    # Compare each nonBCG genome against BCG
    for fasta in glob.glob('{}/*.fasta'.format(bacillus_fasta)):
        accession = splitext(basename(fasta))[0]
        mash_dist = '{0}/{1}.dist'.format(outdir, accession)
        if not os.path.exists(mash_dist):
            run_command([
                'mash', 'dist', '-p', str(args.cpu),
                '{0}/bcg.msh'.format(outdir), fasta
            ], stdout=mash_dist)

        distances = []
        with open(mash_dist, 'r') as mash_handle:
            for line in mash_handle:
                distances.append(float(line.split('\t')[2]))

        if max(distances) <= max(within_bcg_distances):
            # Genome should be included in the group
            print("{0} should be in the group ({1} <= {2})".format(
                genomes[accession]['description'], max(distances),
                max(within_bcg_distances)
            ), file=sys.stderr)
            genomes[accession]['is_bcg'] = 'True'
            genomes[accession]['description'] = (
                '{0}, novel BCG member ({1} <= {2})'
            ).format(genomes[accession]['description'], max(distances),
                     max(within_bcg_distances))

        distances.clear()

    # Output the updated summary
    print('#{0}'.format('\t'.join(cols)))
    for accession, vals in genomes.items():
        print('{0}'.format(
            '\t'.join([vals[c] for c in cols])
        ))

    run_command(['rm', '-rf', bcg_fasta])
    run_command(['rm', '-rf', bacillus_fasta])

