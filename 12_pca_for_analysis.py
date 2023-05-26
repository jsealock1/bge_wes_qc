# PCA on passing scz case/control samples for analysis
import hail as hl 

MT="gs://schema_jsealock/pumas/twist_intervals/pumas_bge_wave_1_qc_out_scz_and_controls.mt"
ANCESTRY = "gs://schema_jsealock/pumas/twist_intervals/pca/bge_wave1_ancestry_assignments_min_prob_0.7.txt"

mt = hl.read_matrix_table(MT)
ancestry = hl.import_table(ANCESTRY, key="Sample")

# filter meta and ancestry to samples to keep
amr = ancestry.filter(ancestry.Population=="amr")
amr_neuromex_paisa = amr.filter((amr.pop=="PUMAS_PAISA_Psychosis_BGE") | (amr.pop=="NeuroMex"))

# filter to amr & neuromex or paisa
mt = mt.filter_cols(hl.is_defined(amr_neuromex_paisa[mt.col_key]))

mt.count_cols()

### PCA
PCA_SNPS = "gs://schema_jsealock/purcell5k_grch38_liftover_2021-09-14.interval_list"
PCA_OUT = "gs://schema_jsealock/pumas/twist_intervals/pca/analysis_pca/neuromex_paisa_amr"

pca_snps = hl.import_locus_intervals(PCA_SNPS, reference_genome="GRCh38")
mt = mt.filter_rows(hl.is_defined(pca_snps[mt.locus]), keep=True)
mt = mt.select_rows('rsid')
mt = mt.select_entries('DP', 'GQ', 'GT')

eigenvalues, score_table, loading_table = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=True)
score_table.flatten().export(f'{PCA_OUT}_pc_scores.txt')
