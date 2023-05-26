## sex_check and ibd filtering
# hailctl dataproc start jsealock-schema --project daly-neale-sczmeta --num-workers 20 --num-preemptible-workers 100 --max-idle 15m 
# hailctl dataproc submit jsealock-schema /Users/juliasealock/Desktop/schema/pumas/07_sex_check.py 


MT = 'gs://schema_jsealock/pumas/twist_intervals/pumas_prefiltered_variants_lcrs_vsqr_intervals.mt'
PRUNED_CHRX_VARIANTS = "gs://schema_jsealock/pumas/sex_check/bge_pumas_pruned_plink_chrX.prune.in"

INITIAL_SAMPLES = "gs://schema_jsealock/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_passing_samples_twist.txt"
IMPUTESEX_TABLE = 'gs://schema_jsealock/pumas/twist_intervals/sex_check/imputesex.ht'
IMPUTESEX_FILE = 'gs://schema_jsealock/pumas/twist_intervals/sex_check/imputesex.tsv'
Y_NCALLED = 'gs://schema_jsealock/pumas/twist_intervals/sex_check/ycalled.tsv'

import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/")

## impute sex
print("load data")
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, key='s')

mt = hl.read_matrix_table(MT)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))

# n = mt.count()

# print('n samples:')
# print(n[1])
# print('n variants:')
# print(n[0])

imputed_sex_default = hl.impute_sex(mt.GT) # is_female: {False: 7819, True: 8897, None: 3142}
imputed_sex_custom = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6) # is_female {False: 7828, True: 12030}

# pruned chrx variants
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
ht_pruned_chrx_variants = ht_pruned_chrx_variants.annotate(**hl.parse_variant(ht_pruned_chrx_variants.f0, reference_genome='GRCh38'))
ht_pruned_chrx_variants = ht_pruned_chrx_variants.key_by(ht_pruned_chrx_variants.locus, ht_pruned_chrx_variants.alleles)
mt2 = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))

imputed_sex_default = hl.impute_sex(mt2.GT) # {False: 7825, True: 10291, None: 1742}
imputed_sex_custom = hl.impute_sex(mt2.GT, female_threshold=0.6, male_threshold=0.6) # {False: 7834, True: 12024}

## use custom imputed sex from pruned variants
mt = mt.annotate_cols(impute_sex = imputed_sex_custom[mt.s])

mt.cols().select('impute_sex').flatten().export(IMPUTESEX_FILE)
# Want to change this to reflect the dataset that I have.
mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

# Determine non-missing allele count on the y.
mt = hl.read_matrix_table(MT)
mt = mt.select_entries(mt.GT).repartition(512)

mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='y_qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.y_qc.n_called).export(Y_NCALLED)


