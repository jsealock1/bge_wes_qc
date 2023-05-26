## filter on sample qc
library(tidyr)
library(dplyr)

# gsutil cp gs://schema_jsealock/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_all_variants_2023-04-04.tsv.bgz .

SAMPLE_QC_PATH = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_all_variants_2023-04-04.tsv.bgz"
PCA_FILTERED_SAMPLES = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/pca/bge_wave1_ancestry_assignments_min_prob_0.7.txt"
COHORT_NAME = "BGE_Wave1"
SAMPLE_METRICS = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/sample_qc/Callset1_PAISA_NeuroMex_QIMR_Kenya_METRICS.csv"
OUT_DIR = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/sample_qc/"

setwd(OUT_DIR)

# subset to samples passing ancestry filter 
dat = read.delim(SAMPLE_QC_PATH)
samples = read.delim(PCA_FILTERED_SAMPLES)
colnames(samples)[4] <- "cohort"
dat$cohort <- NULL
dat = merge(dat, samples, by.x="s", by.y="Sample")
dat$Probability <- NULL

# summarize sample qc metrics
colnames(dat) <- gsub("sample_qc.", "", colnames(dat))

# merge with contam and chimeric percentages 
metrics = read.csv(SAMPLE_METRICS, header=T)
metrics = metrics[c(6,30,32)]
colnames(metrics) <- c("s", "chimeras", "contamination")
dat = merge(metrics, dat, by="s")

## filtering to remove samples with 4 mads away from median within cohort and pop
dat$cohort_pop = paste0(dat$cohort, "_", dat$Population)
fields <- c("n_insertion", "n_deletion", "n_snp", "r_het_hom_var", "r_ti_tv", "r_insertion_deletion")
cohorts <- unique(dat$cohort_pop)

MAD_THRESH = 4
ids <- NA
out = NULL

for (i in 1:length(cohorts)){
  nperm <- dat[dat$cohort_pop == cohorts[i],]
  n = nrow(nperm)
  for (field in fields) {
    # Compute MAD and median
    MAD <- mad(as.numeric(nperm[, field]))
    MED <- median(as.numeric(nperm[, field]))
    # print(sprintf("MAD: %f | MED: %f", MAD, MED))
    # Filter to only entries within 4 MAD of the median
    keep <- nperm[(abs((as.numeric(nperm[, field]) - MED)) < MAD_THRESH*MAD),]$s
    n_removed = sum(!(abs((as.numeric(nperm[, field]) - MED)) < MAD_THRESH*MAD))
    df = data.frame(cohorts[i], field, MAD, MED, n, n_removed)
    out = rbind(out, df) 
  }
  ids <- unique(c(ids, keep))
}

length(ids)
# [1] 20422

# hard filters
### eventually change call rate filter to n bases with gq > 20 / n called
call_rate_t = 0.85
chimeras_t = 15
contam_t = 2

for (i in 1:length(cohorts)){
  nperm <- dat[dat$cohort_pop == cohorts[i],]
  n = nrow(nperm)
  # Compute MAD and median
  field = "call_rate"
  MAD <- mad(as.numeric(nperm[, field]))
  MED <- median(as.numeric(nperm[, field]))
  keep <- subset(nperm, call_rate >= call_rate_t)$s
  n_removed = nrow(subset(nperm, call_rate < call_rate_t))
  df = data.frame(cohorts[i], field, MAD, MED, n, n_removed)
  out = rbind(out, df) 
  # ids1 <- unique(c(ids, keep))
}

for (i in 1:length(cohorts)){
  nperm <- dat[dat$cohort_pop == cohorts[i],]
  n = nrow(nperm)
  # Compute MAD and median
  field = "contamination"
  MAD <- mad(as.numeric(nperm[, field]))
  MED <- median(as.numeric(nperm[, field]))
  keep <- subset(nperm, contamination < contam_t)$s
  n_removed = nrow(subset(nperm, contamination >= contam_t))
  df = data.frame(cohorts[i], field, MAD, MED, n, n_removed)
  out = rbind(out, df) 
  # ids2 <- unique(c(ids, keep))
}

for (i in 1:length(cohorts)){
  nperm <- dat[dat$cohort_pop == cohorts[i],]
  n = nrow(nperm)
  # Compute MAD and median
  field = "chimeras"
  MAD <- mad(as.numeric(nperm[, field]))
  MED <- median(as.numeric(nperm[, field]))
  keep <- subset(nperm, chimeras < chimeras_t)$s
  n_removed = nrow(subset(nperm, chimeras >= chimeras_t))
  df = data.frame(cohorts[i], field, MAD, MED, n, n_removed)
  out = rbind(out,
   df) 
  # ids3 <- unique(c(ids, keep))
}


dat1 = dat[(dat$s %in% ids),]
dat_filt = subset(dat1, call_rate >= call_rate_t & contamination < contam_t & chimeras < chimeras_t)
nrow(dat_filt)
# [1] 19858
write.csv(out, "BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_median_mad_n_filtered_and_hard_filters_twist.csv", row.names=F, quote=F)
dat_filt = dat_filt[c(1,31,32)]
write.table(dat_filt, "BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_passing_samples_twist.txt", col.names=T, row.names=F, quote=F, sep="\t")
