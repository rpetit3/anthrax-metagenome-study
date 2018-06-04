#!/usr/bin/env nextflow
params.help = null
params.outdir = null
params.fasta = null
params.clear_cache_on_success = true
params.cpu = 1
if (params.help) {
    print_usage()
    exit 0
}

check_input_params()

// Set some global variables
name = params.name
outdir = params.outdir ? params.outdir : './'
cpu = params.cpu

/* ==== BEGIN ANNOTATION ==== */
process annotation {
    cache false
    publishDir outdir, mode: 'copy', overwrite: true
    input:
        file fasta from Channel.value(file(params.fasta))
    output:
        file '*.gff'
    shell:
        '''
        prokka --cpus !{cpu} --outdir ./ --force --prefix !{name} \
               --locustag !{name} --quiet !{fasta}
        '''
}
/* ==== END ANNOTATION ==== */

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
    log.info 'Prokka Annotation Pipeline for Roary'
    log.info ''
    log.info 'Required Options:'
    log.info '    --fasta     FASTA   Input FASTA to annotate'
    log.info '    --name      STR     A name to give the run'
    log.info ''
    log.info 'Optional:'
    log.info '    --outdir     DIR  Directory to write results to. (Default ./${NAME})'
    log.info '    --help            Show this message and exit'
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow annotate-genome.nf --fasta input.fasta --name saureus'
}
def check_input_params() {
    error = false
    if (!params.name) {
        log.info('A name is required to continue. Please use --name')
        error = true
    }

    if (!params.fasta) {
        log.info('A reference FASTA is required. Please use --fasta')
        error = true
    } else if (!file(params.fasta).exists()) {
        log.info('Invailid input (--fasta), please verify "' + params.fasta + '"" exists.')
        error = true
    }

    if (error) {
        log.info('See --help for more information')
        exit 1
    }
}
