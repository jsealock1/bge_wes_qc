## pre filters variant QC 
## run after interval filtering, removing LCRs, and failing VQSR and before PCA & sample_qc

## input = matrix table filtered for intervals, removing lcrs, and failing vqsr 

# hailctl dataproc start jsealock-schema --project daly-neale-sczmeta --num-workers 20 --num-preemptible-workers 100 --max-idle 15m --worker-machine-type n1-highmem-8
# hailctl dataproc submit jsealock-schema /Users/juliasealock/Desktop/schema/pumas/twist_intervals/01_pumas_prevariant_qc_sample_qc.py 
# hailctl dataproc stop jsealock-schema

import hail as hl

hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

VCF_PATH = 'gs://bge_callset_paisa_qimr_neuromex_kenyapsych/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych/part_one_outputs/chr*/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_chr*.hard_filtered_with_genotypes.vcf.gz'
# gs://neale-bge/bge-wave-0.vds
META_DATA_PATH = "gs://schema_jsealock/pumas/bge_wave1_meta_data_all_with_phenotype.txt"
SAMPLE_QC_OUT = "gs://schema_jsealock/pumas/twist_intervals/sample_qc/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych_sample_qc_all_variants_2023-04-04.tsv.bgz"

TARGET_INTERVALS = "gs://schema_jsealock/pumas/Twist_Alliance_Clinical_Research_Exome_Covered_Targets_hg38-34.9MB.bed" ## no padding
LCR_PATH = "gs://schema_jsealock/LCRFromHengHg38.bed"
VQSR_VCF = 'gs://bge_callset_paisa_qimr_neuromex_kenyapsych/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych/part_two_outputs/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych.filtered.*.vcf.gz'
PCA_SNPS = "gs://schema_jsealock/purcell5k_grch38_liftover_2021-09-14.interval_list"
PCA_OUT = "gs://schema_jsealock/pumas/twist_intervals/pca/pumas_pca"

# mt = hl.read_matrix_table(MT_PATH)
mt = hl.import_vcf(VCF_PATH, reference_genome='GRCh38', force_bgz=True)

mt = mt.filter_rows(mt.alleles.length() <= 6)
mt = hl.split_multi_hts(mt)



print("pre filter variants")
mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
                # ((mt.AD[0] / mt.DP) < 0.8) | # Has to be removed because allele depth no longer defined for hom ref calls.
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | 
                ((mt.AD[1] / mt.DP) < 0.2) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)


intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr1:START-chr22:END', 'chrX:START-chrX:END', 'chrY:START-chrY:END']]
mt = hl.filter_intervals(mt, intervals)


print("vqsr filter")
part2mt = hl.import_vcf(VQSR_VCF, reference_genome='GRCh38', force_bgz=True)
part2mt = part2mt.annotate_rows(fail_VQSR = hl.len(part2mt.filters) != 0)
passed = part2mt.filter_rows(part2mt.fail_VQSR, keep=False).rows()
mt = mt.filter_rows(hl.is_defined(passed[mt.row_key]), keep=True)

#filter out LCRs 
print("lcr filter")
lcr_intervals = hl.import_locus_intervals(LCR_PATH, reference_genome='GRCh38', skip_invalid_intervals=True)
mt = mt.filter_rows(hl.is_defined(lcr_intervals[mt.locus]), keep=False)

# interval filter
print("interval filter")
intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome="GRCh38")
mt = mt.filter_rows(hl.is_defined(intervals[mt.locus]), keep=True)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

# make into minimal representation 
print("min rep")
mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles))
mt = mt.key_rows_by('min_rep')
mt = mt.drop('locus', 'alleles')

mt = mt.annotate_rows(locus = mt.min_rep.locus,
                      alleles = mt.min_rep.alleles)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.drop('min_rep')

print("n variants after qc")
n = mt.count_rows()
print(n)
# 2205163

print("write prefiltered mt")
mt.write('gs://schema_jsealock/pumas/twist_intervals/pumas_prefiltered_variants_lcrs_vsqr_intervals.mt', overwrite=True)


print("run sample qc")
sample_qc = hl.sample_qc(mt)

print("save sample qc")
sample_qc.cols().select('sample_qc').flatten().export(output=SAMPLE_QC_OUT)



print("run pca")

## pca
print("pca set up")
pca_snps = hl.import_locus_intervals(PCA_SNPS, reference_genome="GRCh38")

##load ref data
gnomad1 = hl.experimental.load_dataset(name='gnomad_hgdp_1kg_subset_dense',
                                  version='3.1.2',
                                  reference_genome='GRCh38',
                                  region='us',
                                  cloud='gcp') # https://hail.is/docs/0.2/experimental/index.html#hail.experimental.load_dataset
gnomad = gnomad1.annotate_cols(Dataset = 'gnomad_hgdp_1kg')
# Take subset of entry fields 
gnomad = gnomad.select_entries('DP', 'GQ', 'GT')
# Take only sites passing VQSR
gnomad = gnomad.filter_rows(gnomad.filters == hl.empty_set(hl.tstr), keep = True)
# Drop other column/global annotations
gnomad = gnomad.select_cols(gnomad.Dataset)
gnomad = gnomad.select_rows('rsid')
gnomad = gnomad.select_globals()

# take purcell 5k snps from both datasets 
gnomad = gnomad.filter_rows(hl.is_defined(pca_snps[gnomad.locus]), keep=True)
mt = mt.filter_rows(hl.is_defined(pca_snps[mt.locus]), keep=True)

# take subset of pumas data to merge with ref data
mt = mt.annotate_cols(Dataset = 'PUMAS')
mt = mt.select_rows('rsid')
mt = mt.select_entries('DP', 'GQ', 'GT')

# combine datasets
mt = mt.union_cols(gnomad)

# Read and repartition
print("write temp1")
mt.write(f'{PCA_OUT}_PCA_temp1.mt', overwrite=True)
print("temp 1 read")
mt = hl.read_matrix_table(f'{PCA_OUT}_PCA_temp1.mt')
mt = mt.repartition(50)
print("temp 2 write") 
mt.write(f'{PCA_OUT}_PCA_temp2.mt', overwrite=True)
print("temp 2 read") 
mt = hl.read_matrix_table(f'{PCA_OUT}_PCA_temp2.mt')

print("save samples lists")
mt.cols().rename({'s' : 'Sample'}).export(f'{PCA_OUT}_PCA_samples.txt') 

pop = gnomad1.select_cols(gnomad_population_inference.pop, hgdp_tgp_meta.project).cols()
pop = pop.select_globals()
pop = pop.filter(pop.project=='synthetic_diploid_truth_sample', keep=False)
pop.flatten().export(f'{PCA_OUT}_gnomad_hgdp_1kgp_pop_labels.txt')


# print("run pca")
# # run pca
eigenvalues, score_table, loading_table = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
score_table.flatten().export(f'{PCA_OUT}_pc_scores.txt')



