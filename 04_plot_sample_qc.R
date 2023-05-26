## plot sample qc metrics

SAMPLE_QC_PATH = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_all_variants_2023-04-04.tsv.bgz"
SAMPLE_QC_PASSING_SAMPLES = 'BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_passing_samples_twist.txt'
PCA_FILTERED_SAMPLES = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/pca/bge_wave1_ancestry_assignments_min_prob_0.7.txt"
COHORT_NAME = "BGE_Wave1"
SAMPLE_METRICS = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/sample_qc/Callset1_PAISA_NeuroMex_QIMR_Kenya_METRICS.csv"
OUT_DIR = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/sample_qc"

setwd(OUT_DIR)

# subset to samples passing ancestry filter 
dat = read.delim(SAMPLE_QC_PATH)
samples = read.delim(PCA_FILTERED_SAMPLES)
passing = read.delim(SAMPLE_QC_PASSING_SAMPLES, header=T)

colnames(samples)[4] <- "cohort"
dat$cohort <- NULL
dat = merge(dat, samples, by.x="s", by.y="Sample")
dat$Probability <- NULL
dat = dat[(dat$s %in% passing$s),]

# summarize sample qc metrics
colnames(dat) <- gsub("sample_qc.", "", colnames(dat))

# merge with contam and chimeric percentages 
metrics = read.csv(SAMPLE_METRICS, header=T)
metrics = metrics[c(6,30,32)]
colnames(metrics) <- c("s", "chimeras", "contamination")
dat = merge(metrics, dat, by="s")
dat$cohort = ifelse(dat$cohort == 'KenyaPsychosis_Koenen_Mamah_BGE', 'Kenya Psychosis',
                ifelse(dat$cohort == 'QIMR_Medland_Bipolar_BGE', 'QIMR Cases',
                  ifelse(dat$cohort == 'NeuroMex', 'NeuroMex',
                    ifelse(dat$cohort == 'QIMR_Whiteman_Controls_BGE', 'QIMR Controls',
                      ifelse(dat$cohort == 'PUMAS_PAISA_Psychosis_BGE', 'Paisa', 'NA')))))
dat$cohort = factor(dat$cohort, levels=c("QIMR Controls", "QIMR Cases", "Paisa", "NeuroMex", "Kenya Psychosis"))
dat$Population = toupper(dat$Population)
dat$Population = factor(dat$Population, levels=rev(c("FIN","MID","SAS",'NFE',"AFR","AMR")))
dat = dat[order(rev(dat$Population)),]
boxplot_fun = function(dat=dat, aes=aes, label=label){
	p = ggplot(dat, aes) +
	    geom_boxplot(outlier.shape = NA) + 
	    geom_jitter(width=0.2, height=0, alpha=0.5, aes(color=Population)) + 
	    coord_flip() +
	    theme_classic() +
	    theme(axis.title.x = element_text(margin = ggplot2::margin(t=10))) +
	    # guides(color="none") +
	    ylab(paste0(label)) +
	    xlab("Cohort") +
	    theme(axis.text=element_text(size=10),
	        axis.title=element_text(size=10,face="bold")) +
	    scale_color_d3() 
return(p)
}

chimeras = boxplot_fun(dat=dat, aes(x=cohort, y=chimeras), label="% Chimeric Reads")
contamination = boxplot_fun(dat=dat, aes(x=cohort, y=contamination), label="% Contamination")
call_rate = boxplot_fun(dat=dat, aes(x=cohort, y=call_rate), label="Call Rate")

titv = boxplot_fun(dat=dat, aes(x=cohort, y=r_ti_tv), label="Ti/Tv Ratio")
del = boxplot_fun(dat=dat, aes(x=cohort, y=n_deletion), label="N Deletion")
ins = boxplot_fun(dat=dat, aes(x=cohort, y=n_insertion), label="N Insertion")
snp = boxplot_fun(dat=dat, aes(x=cohort, y=n_snp), label="N SNPs")
ins_del = boxplot_fun(dat=dat, aes(x=cohort, y=r_insertion_deletion), label="Insertion/Deletion Ratio")
het_hom = boxplot_fun(dat=dat, aes(x=cohort, y=r_het_hom_var), label="Het/Hom Var Ratio")
mean_dp = boxplot_fun(dat=dat, aes(x=cohort, y=dp_stats.mean), label="Mean Depth")
mean_gq = boxplot_fun(dat=dat, aes(x=cohort, y=gq_stats.mean), label="Mean Genotype Quality")

pdf(("sample_qc_pre_filtering_chimeras_contam_call_rate_boxplots_twist_intervals.pdf"), height=10, width=11)
ggarrange(chimeras, contamination, call_rate, ncol=1, common.legend = TRUE, legend="right")
dev.off()

pdf(("sample_qc_pre_filtering_ratios_histograms_twist_intervals.pdf"), height=10, width=11)
ggarrange(titv, ins_del, het_hom, ncol=1, common.legend = TRUE, legend="right")
dev.off()

pdf(("sample_qc_pre_filtering_del_ins_snp_histograms_twist_intervals.pdf"), height=10, width=11)
ggarrange(del, ins, snp, ncol=1, common.legend = TRUE, legend="right")
dev.off()

pdf(("sample_qc_pre_filtering_mean_dp_and_gq_histograms_twist_intervals.pdf"), height=10, width=11)
ggarrange(mean_gq, mean_dp, ncol=1, common.legend = TRUE, legend="right")
dev.off()


