#!/usr/bin/env nextflow
BA_KMERS = "/opt/results/ba-specific-kmers.fasta"
BCG_KMERS = "/opt/results/bcg-specific-kmers.fasta"
LEF_KMERS = "/opt/results/lef-specific-kmers.fasta"
COVERAGES = "/opt/data/coverages.txt"
params.help = null
params.outdir = null
params.fasta = null
params.coverage = COVERAGES
params.replicate = 1
params.clear_cache_on_success = true

if (params.help) {
    print_usage()
    exit 0
}

check_input_params()

// Set some global variables
coverage = null
if (params.coverage.toString().isNumber()) {
    coverage = params.coverage
} else {
    String[] coverages = new File(params.coverage) as String[]
}
name = params.name
outdir = params.outdir ? params.outdir : './'

/* ==== START SIMULATION ==== */

process count_kmers {
    cache false
    publishDir outdir, mode: 'copy', overwrite: true
    input:
        val replicate from Channel.from( 1..params.replicate )
        file fasta from Channel.value(file(params.fasta))
    output:
        file '*.txt.gz'
    shell:
        random_seed = (coverage.toFloat() * replicate.toInteger() * 100).round()
        prefix = name + "-" + coverage + "-" + replicate
        '''
        # Simulate Reads
        art_illumina -l 100 -f !{coverage} -na -ss HS20 -rs !{random_seed} \
                     -i !{fasta} -o !{prefix}

        # Count kmers
        jellyfish count -C -m 31 -s 1M -o !{prefix}.jf !{prefix}.fq

        # Identify group specific kmers
        jellyfish query -s !{BA_KMERS} -o !{prefix}-ba.txt !{prefix}.jf
        gzip --fast !{prefix}-ba.txt
        jellyfish query -s !{BCG_KMERS} -o !{prefix}-bcg.txt !{prefix}.jf
        gzip --fast !{prefix}-bcg.txt
        jellyfish query -s !{LEF_KMERS} -o !{prefix}-lef.txt !{prefix}.jf
        gzip --fast !{prefix}-lef.txt
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
    log.info '    --fasta     FASTA   Input reference FASTA to simulate'
    log.info '    --name      STR     A name to give the run'
    log.info ''
    log.info 'Optional:'
    log.info '    --coverage   FLOAT  Coverage to simulate (Default Pre-Defined Set)'
    log.info '    --replicate  INT    Replicate number for random seed generation (Default 1)'
    log.info '    --outdir     DIR    Directory to write results to. (Default ./${NAME})'
    log.info '    --help              Show this message and exit'
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow illumina-simulations.nf --fasta input.fasta --name saureus'
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
