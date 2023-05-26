setwd("/Users/juliasealock/Desktop/schema/pumas/sex_check")

imputed_sex = "gs://schema_jsealock/pumas/twist_intervals/sex_check/imputesex.tsv"
y_called = "gs://schema_jsealock/pumas/twist_intervals/sex_check/ycalled.tsv"
meta = "/Users/juliasealock/Desktop/schema/pumas/bge_wave1_meta_data_all_with_phenotype.txt"

system(paste0("gsutil cp ", imputed_sex, " ."))
system(paste0("gsutil cp ", y_called, " ."))

imputesex = read.table("imputesex.tsv", header=T)
y_called = read.table('ycalled.tsv', header=T)
meta = read.delim(meta, header=T)

meta_sex = meta[c(3,20)]
impute_sex = imputesex[c(1,2)]

y_called$ycalled = ifelse(y_called$n_called > 0, "TRUE","FALSE")

dat = merge(impute_sex, y_called, by="s")

dat = merge(dat, meta, by.x="s", by.y="SAMPLE_ALIAS")

# dat$bio_sex = ifelse(dat$impute_sex.is_female=="true" & dat$ycalled=="FALSE", "female",
#                 ifelse(dat$impute_sex.is_female=="false" & dat$ycalled=="TRUE", "male", "discordant"))
# table(dat$bio_sex)
# # discordant     female       male 
# #        819      11221       7818
       
# table(dat$bio_sex, dat$SEX)         
#   #            Female Male Unknown
#   # discordant    603   11      52
#   # female       9581   15    1139
#   # male           22 6024     695

# dat$discordant = ifelse(dat$bio_sex=="female" & dat$SEX=="Female", "false",
#                   ifelse(dat$bio_sex=='male' & dat$SEX=="Male", "false", "true"))



discordant1 = subset(dat, (SEX=="Female" & impute_sex.is_female=="false") | (SEX=="Male" & impute_sex.is_female=="true"))
table(discordant1$impute_sex.is_female, discordant1$SEX)
       
  #       Female Male
  # false     34    0
  # true       0   26

discordant1$fail_sex_check = TRUE
discordant1 = discordant1[c(1,4)]
write.table(discordant1, "bge_wave_1_failed_sex_check.tsv", col.names=T, sep="\t", row.names=F, quote=F)

passed_samples = dat[!(dat$s %in% discordant1$s),][c(1,2)]
write.table(passed_samples, "bge_wave_1_passed_sex_check.tsv", col.names=T, sep="\t", row.names=F, quote=F)

meta_failed = meta[(meta$SAMPLE_ALIAS %in% discordant1$s),]
table(meta_failed$COHORT,  meta_failed$PHENOTYPE)
                            
  #                            BP CONTROL MDD SCZ
  # NeuroMex                    4       5   0   7
  # PUMAS_PAISA_Psychosis_BGE  12       7   5   4
  # QIMR_Medland_Bipolar_BGE    2       0   0   0
  # QIMR_Whiteman_Controls_BGE  0      14   0   0

  #################
  ## plots
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

# File locations and plotting locations defined in r_options_BipEx.r
source("r_functions_and_parameters/r_options_BipEx.r")
source("r_functions_and_parameters/pretty_plotting.r")

setwd("/Users/juliasealock/Desktop/schema/pumas/sex_check")

IMPUTESEX_FILE <- "imputesex.tsv"
Y_NCALLED_FILE <- "ycalled.tsv"

SEXCHECK_LIST <- 'sexcheck.remove.sample_list'

df <- fread(IMPUTESEX_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(imputed_sex=as.factor(ifelse(impute_sex.is_female == TRUE, 'Female', 'Male')))

df_y <- fread(Y_NCALLED_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)

df <- merge(df, df_y, by='s')

colors <- pal_d3('category20')(20)[c(1,2)]
fills <- pal_d3('category20')(20)[c(11,12)]

T_impute_sex = 0.6

create_pretty_cumulative(df, aes(impute_sex.f_stat), 'F-statistic', T_impute_sex,
    xlim=c(-T_impute_sex,1.1), title='Cumulative Distribution of F-statistic', save_figure=TRUE, file='05_F_stat_cdf')

pdf("imputed_sex_f_stat_histogram.pdf")
p <- ggplot(df, aes(x=impute_sex.f_stat, fill=imputed_sex)) +
  geom_histogram(binwidth=0.025, alpha=0.8, color='#7f7f7f') +
  scale_fill_manual(values=fills, limits=c('Male', 'Female')) +
  labs(x='X chromosome F-statistic',
       y='Count',
       title='',
       fill='Imputed Sex') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=ggplot2::margin(t=10)),
        axis.title.y=element_text(margin=ggplot2::margin(r=10)),
        plot.title=element_text(hjust=0.5)) +
  geom_vline(xintercept=T_impute_sex, linetype='dashed')
print(p)
dev.off()



## fstat vs reported gender

meta = read.delim("../bge_wave1_meta_data_all_with_phenotype.txt", header=T)
meta_sex = meta[c(3,20)]
colnames(meta_sex) <- c("s", "phenotype.GENDER")
df = merge(df, meta_sex, by="s")

df$phenotype.GENDER <- ifelse(is.na(df$phenotype.GENDER), "Unknown", df$phenotype.GENDER)

## x cohort
p <- ggplot(df, aes(x=impute_sex.f_stat, y=phenotype.LOCATION, colour=phenotype.GENDER)) +
  geom_jitter(width=0, height=0.2, size=1, alpha=0.2, stroke=0.05) + 
  theme_minimal() +
  geom_vline(xintercept=T_impute_sex, linetype='dashed') +
  labs(x='X chromosome F-statistic',
       y='Location',
       color='Reported Sex') 

print(p)
ggsave(paste0(PLOTS, '05_imputesex_scatter_box', '.pdf'), p, width=160, height=90, units='mm')

df_false <- df %>% filter((impute_sex.f_stat > T_impute_sex & phenotype.GENDER == 'Female') | (impute_sex.f_stat < T_impute_sex & phenotype.GENDER == 'Male'))
df_false_plot <- df %>% filter(phenotype.GENDER == 'Unknown' | (impute_sex.f_stat > T_impute_sex & phenotype.GENDER == 'Female') | (impute_sex.f_stat < T_impute_sex & phenotype.GENDER == 'Male'))

# Plots of gender estimates using Giulios plotting method.
pdf("x_f_stat_vs_y_called.pdf")
p <- ggplot(df, aes(x=impute_sex.f_stat, y=impute_sex.n_called, colour=phenotype.GENDER)) +
geom_point(size=0.5) + 
labs(x='X chromosome F-statistic', y='Number of calls in Y', color='Reported Sex') +
scale_color_d3('category10') +
scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
geom_point(data=df_false_plot, aes(x=impute_sex.f_stat, y=impute_sex.n_called), size=0.5) + 
theme_minimal()
print(p)
dev.off()


