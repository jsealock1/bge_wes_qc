## load pca scores and pop/sample labels

library(tidyr)
library(tidyverse)
library(randomForest)
library(ggsci)
library(RColorBrewer)
library(colorRamps)

OUT_DIR = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/pca"
SCORES_PATH = "gs://schema_jsealock/pumas/twist_intervals/pca/pumas_pca_pc_scores.txt"
REF_LABELS = "gs://schema_jsealock/pumas/twist_intervals/pca/gnomad_hgdp_1kgp_pop_labels.txt"
SAMPLES_PATH = "gs://schema_jsealock/pumas/twist_intervals/pca/pumas_pca_PCA_samples.txt"
META = "/Users/juliasealock/Desktop/schema/pumas/bge_wave1_meta_data_all_with_phenotype.txt"

setwd(OUT_DIR)

system(paste0("gsutil cp ", SCORES_PATH, " ."))
# system(paste0("gsutil cp ", REF_LABELS, " ."))
# system(paste0("gsutil cp ", SAMPLES_PATH, " ."))

# load scores 
scores = read.table("pumas_pca_pc_scores.txt", header=T)
scores$scores = gsub("\\[|\\]", "", scores$scores)
scores = separate(data = scores, col = scores, into = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20"), sep = ",")

# load pop labels for ref
gnomad_labels = read.delim("gnomad_hgdp_1kgp_pop_labels.txt", header=T)

# load samples data and merge with meta 
samples = read.table("pumas_pca_PCA_samples.txt", header=T)
pumas = subset(samples, Dataset=="PUMAS")
meta = read.delim(META, header=T)
cohort = meta[c(3,12)]

# create sample, pop, dataset labels for pumas + ref
cohort$project = "BGE Wave 1"
colnames(cohort) <- c("s","pop","project")
labels = rbind(gnomad_labels, cohort)
labels = labels[!is.na(labels$pop),]
# add labels to scores
scores = merge(scores, labels, by="s")
scores = subset(scores, pop!="oth")
scores[,2:21] <- sapply(scores[,2:21],as.numeric)

## RF 
training_data = scores %>% 
  filter(project %in% c('1000 Genomes', 'HGDP')) %>%
  select(pop, PC1:PC10) #  4,040

prediction_data = scores %>% 
  filter(project == 'BGE Wave 1') %>%
  select(s, PC1:PC10)     #  21,107


# RF function to predict ancestry using PCs:
set.seed(42)

# GATK
forest = randomForest(as.factor(pop) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                         data = subset(training_data, 1 == 1),
                         importance = TRUE,
                         ntree = 10000)

# Do the prediction
fit_data = data.frame(predict(forest, prediction_data, type='prob'), Sample = prediction_data$s)

# Do some re-organization voodoo
fit_data = fit_data %>%
  gather(Population, Probability, -Sample) %>%
  group_by(Sample) %>%
  slice(which.max(Probability))

## subset to ancestry predictions with prob > 0.7
fit_cohort = merge(fit_data, cohort, by.x="Sample", by.y="s")
fit_cohort7 = subset(fit_cohort, Probability >= 0.7)
ancestryxcohort = as.data.frame(table(fit_cohort7$Population, fit_cohort7$pop))
ancestryxcohort2 = reshape(ancestryxcohort, idvar="Var2", timevar="Var1", direction="wide")

write.table(fit_cohort, "bge_wave1_ancestry_assignments_all_samples.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(fit_cohort7, "bge_wave1_ancestry_assignments_min_prob_0.7.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(ancestryxcohort2, "bge_wave1_ancestry_assignments_min_prob_0.7_counts.txt", col.names=T, row.names=F, quote=F, sep="\t")

## find failure per cohort
failed = subset(fit_cohort, Probability < 0.7)
failed = as.data.frame(table(failed$pop))
write.table(failed, "bge_wave1_ancestry_assignments_failed_counts.txt", col.names=T, row.names=F, quote=F, sep="\t")


### pretty plots
library(ggsci)

OUT_DIR = "/Users/juliasealock/Desktop/schema/pumas/twist_intervals/pca"
SCORES_PATH = "pumas_pca_pc_scores.txt"
REF_LABELS = "gnomad_hgdp_1kgp_pop_labels.txt"
SAMPLES_PATH = "pumas_pca_PCA_samples.txt"
META = "/Users/juliasealock/Desktop/schema/pumas/bge_wave1_meta_data_all_with_phenotype.txt"

setwd(OUT_DIR)


# load scores 
scores = read.table("pumas_pca_pc_scores.txt", header=T)
scores$scores = gsub("\\[|\\]", "", scores$scores)
scores = separate(data = scores, col = scores, into = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20"), sep = ",")
scores = scores[1:5]

pheno = meta[c(3,24)]

# load pop labels for ref
gnomad_labels = read.delim("gnomad_hgdp_1kgp_pop_labels.txt", header=T)
all_samples = read.table("bge_wave1_ancestry_assignments_all_samples.txt", header=T, sep="\t")

gnomad_labels$phenotype = gnomad_labels$pop
all_samples = all_samples[c(1:4)]
all_samples = merge(all_samples, pheno, by.x="Sample", by.y="SAMPLE_ALIAS")

colnames(all_samples) <- c("s","pop","probability","project", "phenotype")
gnomad_labels$probability = 1
labels = rbind(all_samples, gnomad_labels)
labels$project = ifelse(labels$project=="1000 Genomes" | labels$project=="HGDP", "1KGP + HGDP", ifelse(labels$project == "KenyaPsychosis_Koenen_Mamah_BGE", "Kenya-Psychosis", ifelse(labels$project=="NeuroMex", "NeuroMex", ifelse(labels$project=="PUMAS_PAISA_Psychosis_BGE", "Paisa", ifelse(labels$project=="QIMR_Medland_Bipolar_BGE" | labels$project=="QIMR_Whiteman_Controls_BGE", "QIMR", "NA")))))

# add labels to scores
scores = merge(scores, labels, by="s")
scores = subset(scores, pop!="oth")
scores[,2:5] <- sapply(scores[,2:5],as.numeric)
scores$pop <- toupper(scores$pop)

scores$project = ifelse(scores$probability < 0.7, "failed", scores$project)
scores$project2 = ifelse(scores$project=="1KGP + HGDP", scores$pop, scores$project)
scores$project2 = factor(scores$project2, levels=c("failed","Kenya-Psychosis", "NeuroMex", "Paisa", "QIMR", "AFR", "AMR", "EAS", "MID", "NFE", "SAS"))


pdf(("pc1_vs_pc1_all_samples_ancestry.pdf"), width=8, height=8)
pc1_pc2 = ggplot(scores, aes(x=PC1, y=PC2, color=project2)) + geom_point(alpha = 0.5) + theme_bw() + 
	labs(color="BGE Project or Reference Ancestry") +
	theme(text = element_text(size=14), legend.text=element_text(size=12))
print(pc1_pc2)
dev.off()


pumas = subset(scores, project!="1KGP + HGDP")
ref = subset(scores, project=="1KGP + HGDP")

pumas$phenotype = factor(pumas$phenotype, levels=c("failed","Kenya-Psychosis","NeuroMex","Paisa","QIMR"))
ref$phenotype = factor(ref$phenotype, levels = c("AFR", "AMR", "EAS", "MID", "NFE", "SAS"))

# pal = c("Kenya-Psychosis" = "#E64B35FF", "NeuroMex" = "#4DBBD5FF", "Paisa" = "00A087FF", "QIMR" = "3C5488FF", "AFR" = "#F39B7FFF", "AMR" = "#8491B4FF", "EAS" = "#91D1C2FF", "MID" = "#DC0000FF", "NFE" = "#7E6148FF", "SAS" = "#B09C85FF")

pdf(("pc1_vs_pc2_all_pca_all_samples_ancestry1.pdf"), width=8, height=8)
ggplot() + 
     geom_point(data=ref, aes(x=PC1, y=PC2, fill=project2, color=project2), alpha=0.3, size=4, pch = 21, stroke=NA)+ 
     geom_point(data=pumas, aes(x=PC1, y=PC2, color=project)) + 
     scale_color_manual("phenotype", values=c("failed" = "red", "Kenya-Psychosis" = "#E64B35FF", "NeuroMex" = "#4DBBD5FF", "Paisa" = "#00A087FF", "QIMR" = "#3C5488FF")) +
     scale_fill_manual("project2", values=c("AFR" = "#F39B7FFF", "AMR" = "#8491B4FF", "EAS" = "#91D1C2FF", "MID" = "#DC0000FF", "NFE" = "#7E6148FF", "SAS" = "#B09C85FF")) + 
     theme_bw() + 
     theme(text = element_text(size=14), legend.text=element_text(size=12)) + guides(fill=guide_legend(title="Ancestry")) + guides(color=guide_legend(title="PUMAS Cohort"))
dev.off()


passing = subset(pumas, project!="failed")
pdf(("pc1_vs_pc2_all_pca_all_passing_samples_ancestry1.pdf"), width=8, height=8)
ggplot() + 
     geom_point(data=ref, aes(x=PC1, y=PC2, fill=project2, color=project2), alpha=0.3, size=4, pch = 21, stroke=NA)+ 
     geom_point(data=passing, aes(x=PC1, y=PC2, color=project)) + 
     scale_color_manual("phenotype", values=c("Kenya-Psychosis" = "#E64B35FF", "NeuroMex" = "#4DBBD5FF", "Paisa" = "#00A087FF", "QIMR" = "#3C5488FF")) +
     scale_fill_manual("project2", values=c("AFR" = "#F39B7FFF", "AMR" = "#8491B4FF", "EAS" = "#91D1C2FF", "MID" = "#DC0000FF", "NFE" = "#7E6148FF", "SAS" = "#B09C85FF")) + 
     theme_bw() + 
     theme(text = element_text(size=14), legend.text=element_text(size=12)) + guides(fill=guide_legend(title="Ancestry")) + guides(color=guide_legend(title="PUMAS Cohort"))
dev.off()



# subset to scz cases/controls 
remove = c("BP", "PSYCHOSIS", "OTHER", "MDD")
scz = scores[!(scores$phenotype %in% remove),]
scz$phenotype = toupper(scz$phenotype)
scz$phenotype = ifelse(scz$phenotype=="SCZ","Schizophrenia", ifelse(scz$phenotype=="CONTROL", "Control", scz$phenotype))
scz = subset(scz, project!="QIMR")

ref = subset(scz, project=="1KGP + HGDP")
scz2 = subset(scz, project=="NeuroMex" | project=="Paisa")

scz2$phenotype = factor(scz2$phenotype, levels=c("Schizophrenia", "Control"))
ref$phenotype = factor(ref$phenotype, levels = c("AFR", "AMR", "EAS", "MID", "NFE", "SAS"))

## pal from npg colors in ggsci
pal = c("Schizophrenia" = "red", "Control" = "darkblue", "AFR" = "#F39B7FFF", "AMR" = "#8491B4FF", "EAS" = "#91D1C2FF", "MID" = "#DC0000FF", "NFE" = "#7E6148FF", "SAS" = "#B09C85FF")


pdf(("pc1_vs_pc1_all_pca_scz_case_control_with_ref_ancestry1.pdf"), width=8, height=8)
ggplot() + 
     geom_point(data=ref, aes(x=PC1, y=PC2, fill=project2, color=project2), alpha=0.3, size=3, pch = 21, stroke=NA)+ 
     geom_point(data=scz2, aes(x=PC1, y=PC2, color=phenotype)) + 
     scale_color_manual("phenotype", values=c("Schizophrenia" = "red", "Control" = 'darkblue')) +
     scale_fill_manual("project2", values=c("AFR" = "#F39B7FFF", "AMR" = "#8491B4FF", "EAS" = "#91D1C2FF", "MID" = "#DC0000FF", "NFE" = "#7E6148FF", "SAS" = "#B09C85FF")) + 
     theme_bw() + 
     theme(text = element_text(size=14), legend.text=element_text(size=12)) + guides(fill=guide_legend(title="Ancestry")) + guides(color=guide_legend(title="Case Status"))
dev.off()


