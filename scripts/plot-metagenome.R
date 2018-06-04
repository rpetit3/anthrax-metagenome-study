library(dplyr)
library(ggplot2)
options(scipen=999)

ba_samples = read.table(
    "/home/rpetit/ba-paper/analysis/public-data/nyc-count-summary-sample-ba.txt",
    header=TRUE, sep="\t"
)

bcg_samples = read.table(
    "/home/rpetit/ba-paper/analysis/public-data/nyc-count-summary-sample-bcg.txt",
    header=TRUE, sep="\t"
)

lef_samples = read.table(
    "/home/rpetit/ba-paper/analysis/public-data/nyc-count-summary-sample-lef.txt",
    header=TRUE, sep="\t"
)
nrow(lef_samples[lef_samples$hit > 0,])

ba_samples = merge(
    ba_samples,
    data.frame(
        run=bcg_samples$run,
        bcg_mean_coverage=bcg_samples$kmer_cov_mean
    ),
    by=c('run')
)

# BA, ba coverage vs bcg coverage
ggplot(ba_samples, aes(bcg_mean_coverage, kmer_cov_mean)) +
    geom_point(aes(color=is_ba)) +
    ggtitle("BCG kmer coverage vs BA kmer coverage")

nrow(ba_samples[ba_samples$hit > 0,])

# BA, ba hits vs bcg coverage
ggplot(ba_samples, aes(bcg_mean_coverage, hit)) +
    geom_point(aes(color=is_ba)) +
    ggtitle("BCG kmer coverage vs BA kmer hits")

# BA, ba hits vs ba coverage
ggplot(ba_samples, aes(kmer_cov_mean, hit)) +
    geom_point(aes(color=is_ba)) +
    ggtitle("BCG kmer coverage vs BA kmer hits")


# BA KMER hits
ba_kmers = read.table(
    "/home/rpetit/ba-paper/analysis/public-data/nyc-count-summary-kmer-ba.txt",
    header=TRUE, sep="\t"
)

table(ba_kmers$hamming_distance)

hd_hits = ba_kmers %>%
    group_by(hamming_distance) %>%
    summarise(hit = sum(hit))

hd_kmer_hit = ba_kmers %>%
    group_by(hamming_distance) %>%
    summarise(hit = sum(hit > 0))

hd_miss = ba_kmers %>%
    group_by(hamming_distance) %>%
    summarise(hit = sum(hit == 0))


# BA total hits vs hamming distance
ggplot(hd_hits, aes(hamming_distance, hit)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=hit), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(hd_hits$hamming_distance), max(hd_hits$hamming_distance), by = 1),1)) +
    ggtitle("total kmer counts per hamming distance")
sum(hd_hits$hit)

# BA total individual kmer hits vs hamming distance
ggplot(hd_kmer_hit, aes(hamming_distance, hit)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=hit), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(hd_hits$hamming_distance), max(hd_hits$hamming_distance), by = 1),1)) +
    ggtitle("total individual kmer counts per hamming distance")

sum(hd_kmer_hit$hit)

# BA total hits vs hamming distance
ggplot(hd_miss, aes(hamming_distance, hit)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=hit), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(hd_hits$hamming_distance), max(hd_hits$hamming_distance), by = 1),1)) +
    ggtitle("total individual kmer missed per hamming distance")
sum(hd_miss$hit)


# BCG KMER hits
bcg_kmers = read.table(
    "/home/rpetit/ba-paper/analysis/public-data/nyc-count-summary-kmer-bcg.txt",
    header=TRUE, sep="\t"
)

table(bcg_kmers$hamming_distance)

hd_hits = bcg_kmers %>%
    group_by(hamming_distance) %>%
    summarise(hit = sum(hit))

hd_kmer_hit = bcg_kmers %>%
    group_by(hamming_distance) %>%
    summarise(hit = sum(hit > 0))

hd_miss = bcg_kmers %>%
    group_by(hamming_distance) %>%
    summarise(hit = sum(hit == 0))

# BA total hits vs hamming distance
ggplot(hd_hits, aes(hamming_distance, hit)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=hit), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(hd_hits$hamming_distance), max(hd_hits$hamming_distance), by = 1),1)) +
    ggtitle("BCG total kmer counts per hamming distance")
sum(hd_hits$hit)

# BA total individual kmer hits vs hamming distance
ggplot(hd_kmer_hit, aes(hamming_distance, hit)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=hit), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(hd_hits$hamming_distance), max(hd_hits$hamming_distance), by = 1),1)) +
    ggtitle("BCG total individual kmer counts per hamming distance")
sum(hd_kmer_hit$hit)

# BA total hits vs hamming distance
ggplot(hd_miss, aes(hamming_distance, hit)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=hit), vjust = -0.5) +
    scale_x_continuous(breaks = round(seq(min(hd_hits$hamming_distance), max(hd_hits$hamming_distance), by = 1),1)) +
    ggtitle("BCG total individual kmer missed per hamming distance")
sum(hd_miss$hit)
