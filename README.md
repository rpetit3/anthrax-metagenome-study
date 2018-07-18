[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1283503.svg)](https://doi.org/10.5281/zenodo.1283503)

Fine-scale differentiation between *Bacillus anthracis* and *Bacillus cereus* group signatures in metagenome shotgun data.

# NCBI Query Results
The results for each NCBI query may not return the same results in the future. To counteract this we recorded the results at the time of the query (April 2018).

### All *Bacillus cereus* Group (BCerG) Reference Genomes
```
txid86661[Organism:exp] AND "complete genome"[Title] AND refseq[filter] 
                        AND 3000000:7000000[Sequence Length]
```

The results of this query are available at [data/bcg-summary.txt](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/data/bcg-summary.txt).

### All *Bacillus* (excluding BCerG) Reference Genomes
```
txid1386[Organism:exp] NOT txid86661[Organism:exp] "complete genome"[Title] 
                       AND 3000000:7000000[Sequence Length] AND refseq[filter]
```

The results of this query are available at [data/bacillus-summary.txt](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/data/bacillus-summary.txt).

### pXO1 Plasmids
`pXO1[Title] AND 140000:200000[Sequence Length]`

The results of this query are available at [data/pXO1-accessions.txt](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/data/pXO1-accessions.txt).

### *Bacillus antrhacis* Sequencing Projects
```
genomic[Source] AND random[Selection] AND txid86661[Organism:exp] 
                AND paired[Layout]) AND wgs[Strategy] AND "Illumina HiSeq
```

The results of this query are available at [data/bcg-run-accessions.txt](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/data/bcg-run-accessions.txt)


# Runtime Parameters
Below is a list of runtime parameters for third-party programs.

### ART Sequencing Simulator `art_illumina`
```
art_illumina -l 100 -f $COVERAGE -na -ss HS20 -rs $RANDOM_SEED \
             -i $FASTA_FILE -o $OUTPUT_PREFIX
```

##### Used in the following scripts
- [scripts/illumina-simulations.nf](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/illumina-simulations.nf)
- [scripts/simulate-single-sample.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/simulate-single-sample.py)
- [scripts/subsample-ba-lod.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-ba-lod.py)

### SRA Tools `fastq_dump`
`fastq-dump -I $SRA_ACCESSION`

##### Used in the following scripts
- [scripts/public-data-counts.nf](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/public-data-counts.nf)

### SRA Tools `fastq_dump`
`fastq-dump --split-spot -I --fasta 0 $SRA_ACCESSION`

##### Used in the following scripts
- [scripts/subsample-lethal-factor-cgc.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-lethal-factor-cgc.sh)
- [scripts/subsample-model-set-cgc.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-model-set-cgc.sh)

### Jellyfish `jellyfish count`
`jellyfish count -C -m 31 -s 1M -o $JELLYFISH_DB $FASTA_FILE`

##### Used in the following scripts
- [scripts/download-data.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/download-data.sh)
- [scripts/illumina-simulations.nf](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/illumina-simulations.nf)
- [scripts/public-data-counts.nf](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/public-data-counts.nf)
- [scripts/simulate-single-sample.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/simulate-single-sample.py)
- [scripts/subsample-ba-lod.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-ba-lod.py)
- [scripts/subsample-lethal-factor-cgc.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-lethal-factor-cgc.sh)
- [scripts/subsample-lethal-factor.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-lethal-factor.py)
- [scripts/subsample-model-set-cgc.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-model-set-cgc.sh)
- [scripts/subsample-model-set.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-model-set.py)
- [scripts/subsample-sequences-cgc.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-sequences-cgc.py)
- [scripts/subsample-sequences.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-sequences.py)
- [scripts/subsample-single-sample.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-single-sample.py)

### Jellyfish `jellyfish dump`
`jellyfish dump -c $JELLYFISH_DB`

##### Used in the following scripts
- [scripts/identify-group-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-group-specific-kmers.py)
- [scripts/identify-lef-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-lef-specific-kmers.py)

### Jellyfish `jellyfish query`
`jellyfish query -s $FASTA_FILE $JELLYFISH_DB`

##### Used in the following scripts
- [scripts/identify-group-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-group-specific-kmers.py)
- [scripts/identify-lef-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-lef-specific-kmers.py)
- [scripts/illumina-simulations.nf](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/illumina-simulations.nf)
- [scripts/public-data-counts.nf](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/public-data-counts.nf)
- [scripts/simulate-single-sample.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/simulate-single-sample.py)
- [scripts/subsample-ba-lod.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-ba-lod.py)
- [scripts/subsample-lethal-factor-cgc.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-lethal-factor-cgc.sh)
- [scripts/subsample-lethal-factor.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-lethal-factor.py)
- [scripts/subsample-model-set-cgc.sh](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-model-set-cgc.sh)
- [scripts/subsample-model-set.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-model-set.py)
- [scripts/subsample-sequences-cgc.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-sequences-cgc.py)
- [scripts/subsample-sequences.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-sequences.py)
- [scripts/subsample-single-sample.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/subsample-single-sample.py)

### Mash `mash sketch`
`mash sketch -k 31 -s 100000 -p $NUM_CPU -o $OUTPUT_PREFIX`

##### Used in the following scripts
- [scripts/identify-group-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-group-specific-kmers.py)
- [scripts/identify-lef-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-lef-specific-kmers.py)
- [scripts/identify-novel-bcg.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-novel-bcg.py)

### Mash `mash dist`
`mash dist -p $NUM_CPU $MASH_SKETCH $FASTA_FILE`

##### Used in the following scripts
- [scripts/identify-group-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-group-specific-kmers.py)
- [scripts/identify-lef-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-lef-specific-kmers.py)
- [scripts/identify-novel-bcg.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-novel-bcg.py)

### 31-mer Hamming Distance `blastn`
```
blastn -max_hsps 1 -max_target_seqs 1 -dust no -word_size 7 -outfmt 15 \
       -query $FASTA_FILE -db $BLAST_DB -evalue 10000 -num_threads $NUM_CPU
```

##### Used in the following scripts
- [scripts/identify-group-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-group-specific-kmers.py)
- [scripts/identify-lef-specific-kmers.py](https://github.com/rpetit3/anthrax-metagenome-study/blob/master/scripts/identify-lef-specific-kmers.py)

### Other or Previous Analysis
The following are from our previous analysis available at [nyc-subway-anthrax-study](https://github.com/Read-Lab-Confederation/nyc-subway-anthrax-study/)

#### BWA `bwa mem`
`bwa mem -t $NUM_CPU $REFERENCE $FASTQ_R1 $FASTQ_R2 > $SAM_FILE`

#### samtools `samtools view | sort`
`samtools view -@ 10 -bS $SAM_FILE | samtools sort -@ 10 - $SAMPLE`

#### samtools `index`
`samtools index -b $BAM_FILE`

#### bedtools `genomeCoverageBed`
`genomeCoverageBed -d -ibam $BAM_FILE | gzip --best - > $COVERAGE`

#### BLAST+ `blastn`
```
blastn -query $FASTA_FILE -db $BLAST_NT -num_threads $NUM_CPU -max_target_seqs 5 \
    -outfmt "6 stitle sseqid qseqid qstart qend sstart send evalue bitscore score \
               length pident nident mismatch positive means gapopen ppos qcovs qcovhsp" 
    > $BLAST_OUTPUT
```

#### `bam2fastq`
`bam2fastq -o $FASTQ --no-unaligned $BAM_FILE`

#### FASTX `fastq_to_fasta`
`cat $FASTQ_FILE | fastq_to_fasta -Q33 -n | gzip --best - > $FASTA_OUTPUT`

#### Kraken
`kraken-build --standard --db $DATABASE`

* Implemented https://github.com/DerrickWood/kraken/pull/57 to auto apply dust masking via `dustmasker`
