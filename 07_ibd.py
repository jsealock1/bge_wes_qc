## ibd
# hailctl dataproc start jsealock-schema --project daly-neale-sczmeta --num-workers 20 --num-preemptible-workers 100 --max-idle 15m 
# hailctl dataproc submit jsealock-schema /Users/juliasealock/Desktop/schema/pumas/twist_intervals/07_ibd.py 


import hail as hl
hl.init()

MT = 'gs://schema_jsealock/pumas/twist_intervals/pumas_prefiltered_variants_lcrs_vsqr_intervals.mt'
# IBD_OUTPUT = 'gs://schema_jsealock/pumas/twist_intervals/relatedness/pc_relate_pairs_min_0.177.tsv'
INITIAL_SAMPLES = 'gs://schema_jsealock/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_passing_samples_twist.txt'
PASS_SEX_CHECK = "gs://schema_jsealock/pumas/twist_intervals/sex_check/bge_wave_1_passed_sex_check.tsv"
PRUNED_VARIANTS = "gs://schema_jsealock/purcell5k_grch38_liftover_2021-09-14.interval_list"
PHENOTYPES_TABLE = 'gs://schema_jsealock/pumas/bge_wave1_meta_data_all_with_phenotype.txt'
# PASSING_SAMPLES = 'gs://schema_jsealock/pumas/twist_intervals/relatedness/samples_passing_ibd_filter.tsv'

print("load variant and sample data")
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, key='s')
ht_sex_checked = hl.import_table(PASS_SEX_CHECK, key='s')
pruned_variants_snps = hl.import_locus_intervals(PRUNED_VARIANTS, reference_genome="GRCh38")

sample_annotations = hl.import_table(PHENOTYPES_TABLE, key="SAMPLE_ALIAS")
sample_annotations = sample_annotations.annotate(scz_status = hl.if_else(sample_annotations.PHENOTYPE == "SCZ", hl.bool("TRUE"), hl.bool("FALSE")))

print("load and filter mt")
mt = hl.read_matrix_table(MT)
mt = mt.select_entries(mt.GT).repartition(512)

mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(hl.is_defined(ht_sex_checked[mt.col_key]))

mt = mt.filter_rows(hl.is_defined(pruned_variants_snps[mt.locus]), keep=True)
mt = mt.annotate_cols(scz_status = sample_annotations[mt.s].scz_status).repartition(128).persist()
# {False: 17479, True: 2379}

mt = mt.annotate_cols(cohort = sample_annotations[mt.s].COHORT)
mt = mt.filter_cols((mt.cohort == "NeuroMex") | (mt.cohort=="PUMAS_PAISA_Psychosis_BGE"))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s].PHENOTYPE) #.repartition(128).persist()
mt = mt.filter_cols((mt.phenotype == "SCZ") | (mt.phenotype=="CONTROL"))
## starting 
#### {'CONTROL': 2298, 'SCZ': 2368} 
#### {'NeuroMex': 1791, 'PUMAS_PAISA_Psychosis_BGE': 2875}

######################
## use pc_relate and max_ind_set in hail to preferentially keep scz cases
# https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set

print("run pc_relate")
dataset = mt
pc_rel = hl.pc_relate(dataset.GT, 0.001, k=5, statistics='kin')
pairs = pc_rel.filter(pc_rel['kin'] > 0.177)
IBD_OUTPUT = 'gs://schema_jsealock/pumas/twist_intervals/relatedness/pc_relate_pairs_min_0.177_scz_and_controls_paisa_neuromex.tsv'
pairs.export(IBD_OUTPUT)

# Starting from the above pairs, prune individuals from a dataset until no close relationships remain, preferring to keep cases over controls:

print("filter for max ind set")
## scz 
samples = dataset.cols()
pairs_with_case = pairs.key_by(
    i=hl.struct(id=pairs.i, is_case=samples[pairs.i].scz_status),
    j=hl.struct(id=pairs.j, is_case=samples[pairs.j].scz_status))
def tie_breaker(l, r):
    return hl.if_else(l.is_case & ~r.is_case, -1,
                      hl.if_else(~l.is_case & r.is_case, 1, 0))
related_samples_to_remove = hl.maximal_independent_set(
   pairs_with_case.i, pairs_with_case.j, False, tie_breaker)
result = dataset.filter_cols(hl.is_defined(
	related_samples_to_remove.key_by(
		s = related_samples_to_remove.node.id.s)[dataset.col_key]), keep=False)

n_remove = result.count_cols()
print('n samples kept:')
print(n_remove)

print("export passing samples")
samples_keep = result.cols()
PASSING_SAMPLES = 'gs://schema_jsealock/pumas/twist_intervals/relatedness/samples_passing_ibd_filter_scz_and_controls_paisa_neuromex.tsv'
samples_keep.export(PASSING_SAMPLES)



#### create max independent sets for other phenotypes 

# ## bip
# samples = dataset.cols()
# pairs_with_case = pairs.key_by(
#     i=hl.struct(id=pairs.i, is_case=samples[pairs.i].bp_status),
#     j=hl.struct(id=pairs.j, is_case=samples[pairs.j].bp_status))
# def tie_breaker(l, r):
#     return hl.if_else(l.is_case & ~r.is_case, -1,
#                       hl.if_else(~l.is_case & r.is_case, 1, 0))
# related_samples_to_remove = hl.maximal_independent_set(
#    pairs_with_case.i, pairs_with_case.j, False, tie_breaker)
# result = dataset.filter_cols(hl.is_defined(
#     related_samples_to_remove.key_by(
#        s = related_samples_to_remove.node.id.s)[dataset.col_key]), keep=False)


# ## mdd
# samples = dataset.cols()
# pairs_with_case = pairs.key_by(
#     i=hl.struct(id=pairs.i, is_case=samples[pairs.i].mdd_status),
#     j=hl.struct(id=pairs.j, is_case=samples[pairs.j].mdd_status))
# def tie_breaker(l, r):
#     return hl.if_else(l.is_case & ~r.is_case, -1,
#                       hl.if_else(~l.is_case & r.is_case, 1, 0))
# related_samples_to_remove = hl.maximal_independent_set(
#    pairs_with_case.i, pairs_with_case.j, False, tie_breaker)
# result = dataset.filter_cols(hl.is_defined(
#     related_samples_to_remove.key_by(
#        s = related_samples_to_remove.node.id.s)[dataset.col_key]), keep=False)


# ## psychosis
# samples = dataset.cols()
# pairs_with_case = pairs.key_by(
#     i=hl.struct(id=pairs.i, is_case=samples[pairs.i].psychosis_status),
#     j=hl.struct(id=pairs.j, is_case=samples[pairs.j].psychosis_status))
# def tie_breaker(l, r):
#     return hl.if_else(l.is_case & ~r.is_case, -1,
#                       hl.if_else(~l.is_case & r.is_case, 1, 0))
# related_samples_to_remove = hl.maximal_independent_set(
#    pairs_with_case.i, pairs_with_case.j, False, tie_breaker)
# result = dataset.filter_cols(hl.is_defined(
#     related_samples_to_remove.key_by(
#        s = related_samples_to_remove.node.id.s)[dataset.col_key]), keep=False)
