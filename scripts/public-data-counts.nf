#!/usr/bin/env nextflow
BA_KMERS = "/opt/results/ba-specific-kmers.fasta"
BCG_KMERS = "/opt/results/bcg-specific-kmers.fasta"
LEF_KMERS = "/opt/results/lef-specific-kmers.fasta"
params.help = null
params.outdir = null
params.accession = null
params.cpu = 4
params.clear_cache_on_success = true

if (params.help) {
    print_usage()
    exit 0
}

check_input_params()

// Set some global variables
accession = params.accession
outdir = params.outdir ? params.outdir : './'
cpu = params.cpu

/* ==== START SIMULATION ==== */

process count_kmers {
    cache false
    publishDir outdir, mode: 'copy', overwrite: true

    output:
        file '*.tar'
    shell:
        '''
        # Download Reads
        fastq-dump -I !{accession}

        # Count kmers
        jellyfish count -C -m 31 -s 1M -t !{cpu} -o !{accession}.jf !{accession}.fastq

        # Identify group specific kmers
        jellyfish query -s !{BA_KMERS} -o !{accession}-ba.txt !{accession}.jf
        gzip --fast !{accession}-ba.txt
        jellyfish query -s !{BCG_KMERS} -o !{accession}-bcg.txt !{accession}.jf
        gzip --fast !{accession}-bcg.txt
        jellyfish query -s !{LEF_KMERS} -o !{accession}-lef.txt !{accession}.jf
        gzip --fast !{accession}-lef.txt
        tar cvf !{accession}.tar *.txt.gz
        '''
}

/* ==== END SIMULATION ==== */

workflow.onComplete {
    if (workflow.success == true && params.clear_cache_on_success) {
        // No need to resume completed run so remove cache.
        file('./work/').deleteDir()
    }
    println """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

// Utility Functions
def print_usage() {
    log.info 'Illumina Simulation Pipeline'
    log.info ''
    log.info 'Required Options:'
    log.info '    --accession ACCESSION   ENA/SRA run accession to identify kmers'
    log.info ''
    log.info 'Optional:'
    log.info '    --outdir     DIR  Directory to write results to. (Default ./${NAME})'
    log.info '    --help            Show this message and exit'
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow illumina-simulations.nf --accession RUN_ACCESSION'
}
def check_input_params() {
    error = false
    if (!params.accession) {
        log.info('A ENA/SRA run accession is required to continue. Please use --accession')
        error = true
    }

    if (error) {
        log.info('See --help for more information')
        exit 1
    }
}
