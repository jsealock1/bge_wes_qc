## filter for passing samples from sample qc, independent samples from ibd, and scz case/control 

# hailctl dataproc start jsealock-schema --project daly-neale-sczmeta --num-workers 20 --num-preemptible-workers 100 --max-idle 15m 
# hailctl dataproc submit jsealock-schema /Users/juliasealock/Desktop/schema/pumas/twist_intervals/08_filter_matrix_table_passing_samples.py
import hail as hl
hl.init()

MT = 'gs://schema_jsealock/pumas/twist_intervals/pumas_prefiltered_variants_lcrs_vsqr_intervals.mt'
INITIAL_SAMPLES = 'gs://schema_jsealock/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_passing_samples_twist.txt'
PASS_SEX_CHECK = "gs://schema_jsealock/pumas/twist_intervals/sex_check/bge_wave_1_passed_sex_check.tsv"
PHENOTYPES_TABLE = 'gs://schema_jsealock/pumas/bge_wave1_meta_data_all_with_phenotype.txt'

MT_OUT = 'gs://schema_jsealock/pumas/twist_intervals/pumas_bge_wave_1_qc_out_all_samples.mt' 

MT_SCZ_OUT = 'gs://schema_jsealock/pumas/twist_intervals/pumas_bge_wave_1_qc_out_scz_and_controls.mt'
SCZ_UNRELATED_SAMPLES = 'gs://schema_jsealock/pumas/twist_intervals/relatedness/samples_passing_ibd_filter_scz_and_controls_paisa_neuromex.tsv'

MT_PAISA = "gs://pumas_qced_data/pumas_bge_wave_1_qc_out_paisa_samples.mt "## change this to shared bucket

print("load data")
mt = hl.read_matrix_table(MT)
scz_unrelated = hl.import_table(SCZ_UNRELATED_SAMPLES, key="s")
initial_samples = hl.import_table(INITIAL_SAMPLES, key='s')
ht_sex_checked = hl.import_table(PASS_SEX_CHECK, key='s')

print("filter mt")
mt = mt.filter_cols(hl.is_defined(initial_samples[mt.col_key]))
mt = mt.filter_cols(hl.is_defined(ht_sex_checked[mt.col_key]))

print("annotate with pheno and cohort")
sample_annotations = hl.import_table(PHENOTYPES_TABLE, key="SAMPLE_ALIAS")
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s].PHENOTYPE) #.repartition(128).persist()
mt = mt.annotate_cols(cohort = sample_annotations[mt.s].COHORT)

## create max independent annotations cols for scz, bip, mdd, and psychosis 
print("write full mt")
mt.write(MT_OUT)

## write out paisa subset to shared bucket 
print("write paisa")
paisa = mt.filter_cols(mt.cohort == 'PUMAS_PAISA_Psychosis_BGE', keep=True)
paisa.write(MT_PAISA)

# print("annotate with phenotype and filter for scz/controls")
# scz = mt.filter_cols((mt.phenotype == "SCZ") | (mt.phenotype=="CONTROL"))
# scz = scz.filter_cols(hl.is_defined(scz_unrelated[scz.col_key]))
# # remove qimr controls bc no cases 
# scz = scz.filter_cols(scz.cohort == "QIMR_Whiteman_Controls_BGE", keep=False)
# # total 4619
# # {'NeuroMex': 1780, 'PUMAS_PAISA_Psychosis_BGE': 2839}
# # {'CONTROL': 2264, 'SCZ': 2355}

# print("write scz/control mt")
# scz.write(MT_SCZ_OUT)







