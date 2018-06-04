library(ggplot2)
options(scipen=999)

read_samples <- function(gzip_file, distance) {
    s = read.table(gzfile(gzip_file), header=TRUE)

    # See Olson et al. 2015 frontiers in Genetics for calculations
    s$accuracy = ifelse(s$tp + s$tn > 0, (s$tn + s$fp) / s$total_kmers, 0) # (tp + tn) / total
    s$specificity = ifelse(s$fp + s$tn > 0, s$tn / (s$tn + s$fp), 0) # tn / (fp + tn)
    s$sensitivity = ifelse(s$tp + s$fn > 0, s$tp / (s$tp + s$fn), 0) # tp / (tp + fn)
    s$precision = ifelse(s$tp + s$fp > 0, s$tp / (s$tp + s$fp), 0) # tp / (tp + fp)
    s$fpr = ifelse(s$fp + s$tn > 0, s$fp / (s$tn + s$fp), 0) # fp / (tn + fp)

    # bacilllus group
    s$bacillus = ifelse(
        s$is_bcg == 'True',
        ifelse(s$is_ba=='True', 'ba', 'bcg'),
        'nonbcg'
    )

    s = merge(
        s,
        data.frame(
            sample=distance$accession,
            distance=distance$distance
        ),
        by=c('sample')
    )

    return(s)
}

ba_kmer_info = read.table(
    "/home/rpetit/ba-paper/results/ba-kmer-info.txt",
    header=TRUE, sep="\t"
)

bcg_kmer_info = read.table(
    "/home/rpetit/ba-paper/results/bcg-kmer-info.txt",
    header=TRUE, sep="\t"
)


ba_samples = read_samples(
    "/home/rpetit/ba-paper/analysis/simulations/count-summary-sample-ba.txt.gz",
    ba_kmer_info
)

bcg_samples = read_samples(
    "/home/rpetit/ba-paper/analysis/simulations/count-summary-sample-bcg.txt.gz",
    ba_kmer_info
)

ba_samples = merge(
    ba_samples,
    data.frame(
        sample=bcg_samples$sample,
        simulated_coverage=bcg_samples$simulated_coverage,
        bcg_mean_coverage=bcg_samples$kmer_cov_mean
    ),
    by=c('sample', 'simulated_coverage')
)




lef_samples = read_samples(
    "/home/rpetit/ba-paper/analysis/simulations/count-summary-sample-lef.txt.gz",
    ba_kmer_info
)

# BA Sensitivity (TPR)
ggplot(ba_samples[ba_samples$is_ba == 'True',], aes(simulated_coverage, sensitivity)) +
    geom_point(aes(color=is_ba)) +
    ggtitle("simulated coverage vs BA kmer Sensitivity")

# BCG Sensitivity (TPR)
ggplot(bcg_samples[bcg_samples$is_bcg == 'True',], aes(simulated_coverage, sensitivity)) +
    geom_point(aes(color=is_bcg)) +
    ggtitle("simulated coverage vs BCG kmer Sensitivity")

# BA, ba coverage vs bcg coverage
ggplot(ba_samples[ba_samples$is_ba == 'True',], aes(bcg_mean_coverage, kmer_cov_mean)) +
    geom_point(aes(color=is_ba)) +
    ggtitle("BCG kmer coverage vs BA kmer coverage")

# BCG BA False Positive Rate
ggplot(ba_samples[ba_samples$bacillus == 'bcg',], aes(bcg_mean_coverage, fpr)) +
    geom_point(aes(color=distance)) +
    ggtitle("BCG kmer coverage vs BA fpr, non BA BCG members")

ggplot(ba_samples[ba_samples$bacillus == 'bcg',], aes(simulated_coverage, fpr)) +
    geom_point(aes(color=distance)) +
    ggtitle("simulated coverage vs BA fpr, non BA BCG members")

# BCG BA False Positive Hits
ggplot(ba_samples[ba_samples$bacillus == 'bcg',], aes(bcg_mean_coverage, fp)) +
    geom_point(aes(color=distance)) +
    ggtitle("BCG kmer coverage vs False positive BA count, non BA BCG members")

ggplot(ba_samples[ba_samples$bacillus == 'bcg',], aes(simulated_coverage, fp)) +
    geom_point(aes(color=distance)) +
    ggtitle("simulated coverage vs False positive BA count, non BA BCG members")


# nonBCG BA False Positive Rate
ggplot(ba_samples[ba_samples$bacillus == 'nonbcg',], aes(simulated_coverage, fpr)) +
    geom_point(aes(color=distance)) +
    ggtitle("simulated coverage vs BA fpr, nonBCG bacillus")

# nonBCG BA False Positive Hits
ggplot(ba_samples[ba_samples$bacillus == 'nonbcg',], aes(simulated_coverage, fp)) +
    geom_point(aes(color=distance)) +
    ggtitle("simulated coverage vs BA false positive hits, nonBCG bacillus")
