#! /usr/bin/env python3
# A wrapper for executing Prokka annotation.
import logging
import os
import subprocess
from Bio import Entrez
Entrez.email = 'robert.petit@emory.edu'
DATABASE = 'nucleotide'
RETTYPE = 'fasta'
RETMODE = 'text'


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
        logging.error('END\n')
        raise RuntimeError(err)
    else:
        logging.info('COMMAND: {0}'.format(cmd))
        logging.info('STDOUT: {0}'.format(out))
        logging.info('STDERR: {0}'.format(err))
        logging.info('END\n')
        return [out, err]


def byte_to_string(string):
    string = string.decode("utf-8") if string else ''
    return string


def run_command(cmd, cwd=os.getcwd(), stdout=False, stderr=False):
    """Execute a single command and return STDOUT and STDERR."""
    stdout, stdout_str = output_handler(stdout)
    stderr, stderr_str = output_handler(stderr, redirect='2>')

    p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, cwd=cwd)

    out, err = p.communicate()
    return onfinish_handler(
        '{0} {1} {2}'.format(' '.join(cmd), stdout_str, stderr_str),
        byte_to_string(out), byte_to_string(err), p.returncode
    )


def generate_nextflow(name, fasta, resume, cpu=1):
    cmd = ['./annotate-genome.nf', '--name', name, '--fasta', fasta,
           '--cpu', str(cpu)]

    if resume:
        cmd.append('-resume')

    return cmd


if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='annotate-genome.py',
        conflict_handler='resolve',
        description=('A wrapper for executing Prokka annotation.'))
    parser.add_argument('accession', metavar="ACCESSION", type=str,
                        help=('NCBI accession to retrieve reference FASTA'))
    parser.add_argument(
        '--cpu', metavar="INT", type=str, default="1",
        help='Number of processors to use (Default 1)'
    )
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Tell nextflow to resume the run.')

    args = parser.parse_args()
    name = args.accession
    outdir = '{0}/{1}'.format(os.getcwd(), args.accession)

    # Setup logs
    logging.basicConfig(filename='{0}-annotation.txt'.format(args.accession),
                        filemode='w', level=logging.INFO)

    # Make directory
    run_command(['mkdir', '-p', args.accession])

    # Download FASTA
    fasta = '{0}/{1}/{1}.fasta'.format(os.getcwd(), args.accession)
    with open(fasta, 'w') as fasta_handle:
        efetch = Entrez.efetch(db=DATABASE, id=args.accession,
                               rettype=RETTYPE, retmode=RETMODE)
        fasta_handle.write(efetch.read())
        efetch.close()

    # Run pipeline
    run_command(['cp', '/usr/local/bin/annotate-genome.nf', outdir])
    nextflow = generate_nextflow(
        args.accession, fasta, args.resume, cpu=args.cpu
    )
    run_command(nextflow, cwd=outdir)

    # Tarball and delete directory.nextflow.log
    run_command(['mv', '{0}/{1}.gff'.format(outdir, args.accession), './'])
    run_command(['rm', '-rf', outdir])
