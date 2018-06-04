#! /usr/bin/env python3
# Simulate illumina sequencing and identify group specific kmers.
import logging
import os
import subprocess
from Bio import Entrez
Entrez.email = 'robert.petit@emory.edu'
COVERAGES = "/opt/data/coverages.txt"
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


def generate_nextflow(name, fasta, replicate, resume, coverage):
    cmd = ['./illumina-simulations.nf', '--name', name, '--fasta', fasta,
           '--replicate', replicate]

    if coverage:
        cmd.append('--coverage')
        cmd.append(coverage)

    if resume:
        cmd.append('-resume')

    return cmd


if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='illumina-simulations.py',
        conflict_handler='resolve',
        description=('A wrapper for executing Illumina Simualtion.'))
    parser.add_argument('accession', metavar="ACCESSION", type=str,
                        help=('NCBI accession to retrieve reference FASTA'))

    parser.add_argument('--coverage', metavar="FLOAT", type=str,
                        help='Text file with coverages to simulate.')
    parser.add_argument(
        '--replicates', metavar="INT", type=int, default=20,
        help='Number of replicates to create (Default 20)'
    )
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Tell nextflow to resume the run.')

    args = parser.parse_args()
    name = args.accession
    coverages = []
    if args.coverage:
        if os.path.exists(args.coverage):
            with open(args.coverage, 'r') as input_handle:
                for line in input_handle:
                    coverages.append(line.rstrip())
        else:
            if float(args.coverage):
                coverages.append(args.coverage)
    else:
        with open(COVERAGES, 'r') as input_handle:
            for line in input_handle:
                coverages.append(line.rstrip())

    # Setup logs
    log_file = '{0}{1}-simulation.txt'.format(
        args.accession, "-{0}".format(args.coverage) if args.coverage else ""
    )
    logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO)

    # Make directory
    run_command(['mkdir', '-p', name])

    # Download FASTA
    fasta = '{0}/{1}/{1}.fasta'.format(os.getcwd(), name)
    with open(fasta, 'w') as fasta_handle:
        efetch = Entrez.efetch(db=DATABASE, id=name,
                               rettype=RETTYPE, retmode=RETMODE)
        fasta_handle.write(efetch.read())
        efetch.close()

    if coverages:
        for coverage in coverages:
            outdir = '{0}/{1}/{2}'.format(os.getcwd(), name, coverage)
            run_command(['mkdir', '-p', outdir])
            # Run pipeline
            run_command(
                ['cp', '/usr/local/bin/illumina-simulations.nf', outdir]
            )
            nextflow = generate_nextflow(
                name, fasta, str(args.replicates), args.resume, coverage
            )
            run_command(nextflow, cwd=outdir)
            run_command(['rm', '-rf', "{0}/.nextflow/".format(outdir)])
            run_command(['rm', '-rf', "{0}/.nextflow.log".format(outdir)])
            run_command(
                ['rm', '-rf', "{0}/illumina-simulations.nf".format(outdir)]
            )

        # Tarball and delete directory.nextflow.log
        run_command(['rm', '-rf', fasta])
        run_command(['tar', '-cvf', '{0}-simulation.tar'.format(name), name])
        run_command(['rm', '-rf', name])
