## find all rare variants
# hailctl dataproc start js --project daly-neale-sczmeta --num-workers 10 --num-preemptible-workers 20 --max-idle 15m 
# hailctl dataproc submit js /Users/juliasealock/Desktop/schema/pumas/twist_intervals/11_agg_variant_per_gene.py
# hailctl dataproc stop js

import hail as hl
hl.init()

from pprint import pprint

# Input
MT="gs://schema_jsealock/pumas/twist_intervals/pumas_bge_wave_1_qc_out_scz_and_controls.mt"
SITES_TABLE = 'gs://schema_jsealock/pumas/twist_intervals/annotations/sites_annot.ht'

RESULTS_OUT = "gs://schema_jsealock/pumas/twist_intervals/results"

############
print("")
print("Read in data")
mt = hl.read_matrix_table(MT)
ht_site = hl.read_table(SITES_TABLE)

print("filter to AMR ancestry")
ancestry = hl.import_table("gs://schema_jsealock/pumas/twist_intervals/pca/bge_wave1_ancestry_assignments_min_prob_0.7.txt", key='Sample')
amr = ancestry.filter(ancestry.Population=="amr")
mt = mt.filter_cols(hl.is_defined(amr[mt.col_key]))

print("Annotate with variant annotations")
mt = mt.annotate_rows(site = ht_site[mt.row_key])

OUT_FILE = f'{RESULTS_OUT}/bge_wave1_amr_counts_by_annotation_ac5.tsv'
MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_all_ac5.tsv'
MPC3_MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_MPC_greater_3_AC5.tsv'
MPC2_TO_3_MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_MPC_2-3_AC5.tsv'
MPC2_MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_MPC_at_least_2_AC5.tsv'
MPC_LESS_2_MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_MPC_less_than_2_AC5.tsv'
SYN_MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_synonymous_variants_AC5.tsv'
PTV_MAC5_TABLE = f'{RESULTS_OUT}/bge_wave1_amr_counts_ptv_ac5.tsv'

####
mt = hl.variant_qc(mt)
mt = mt.annotate_rows(
    ptv = mt.site.consequence_category == "pLoF",
    mis1 = mt.site.consequence_category == "other_missense",
    mis3 = mt.site.consequence_category == "damaging_missense",
    syn = mt.site.consequence_category == "synonymous",
    mpc3 = mt.site.MPC > 3,
    mpc2_to_3 = (mt.site.MPC >=2) & (mt.site.MPC <= 3),
    mpc_less_2 = mt.site.MPC <2,
    mpc_at_least_2 = mt.site.MPC >=2,
    isSingleton = (mt.variant_qc.n_non_ref==1),
    AC5=(mt.variant_qc.AC[1] < 6),
    AC10=(mt.variant_qc.AC[1] < 11),
    )

mt2 = mt.group_rows_by(mt.site.vep.worst_csq_for_variant_canonical.gene_symbol).aggregate(
    n_ac5 = hl.agg.count_where(mt.AC5 & mt.GT.is_non_ref()),
    mpc3_ac5 = hl.agg.count_where(mt.AC5 & mt.mpc3 & mt.GT.is_non_ref()),
    mpc2_to_3_ac5 = hl.agg.count_where(mt.AC5 & mt.mpc2_to_3 & mt.GT.is_non_ref()),
    mpc_less_2_ac5 = hl.agg.count_where(mt.AC5 & mt.mpc_less_2 & mt.GT.is_non_ref()),
    mpc_at_least_2_ac5 = hl.agg.count_where(mt.AC5 & mt.mpc_at_least_2 & mt.GT.is_non_ref()),
    syn_ac5 = hl.agg.count_where(mt.AC5 & mt.syn & mt.GT.is_non_ref()),
    n_ptv_ac5 = hl.agg.count_where(mt.AC5 & mt.ptv & mt.GT.is_non_ref())
    )

mt2.n_ac5.export(MAC5_TABLE)
mt2.mpc3_ac5.export(MPC3_MAC5_TABLE)
mt2.mpc2_to_3_ac5.export(MPC2_TO_3_MAC5_TABLE)
mt2.mpc_less_2_ac5.export(MPC_LESS_2_MAC5_TABLE)
mt2.mpc_at_least_2_ac5.export(MPC2_MAC5_TABLE)
mt2.syn_ac5.export(SYN_MAC5_TABLE)
mt2.n_ptv_ac5.export(PTV_MAC5_TABLE)




print("filter to genes and find total n ac5")

pli = hl.import_table("gs://schema_jsealock/gnomad_subset_qc/annotations/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", key="gene", impute=True)

mt2 = mt2.filter_rows(hl.is_defined(pli[mt2.row_key]))
mt2 = mt2.annotate_cols(total_n_ac5 = hl.agg.sum(mt2.n_ac5))

print("filter to high pli")
pli = pli.filter(pli.pLI > 0.9)
mt2 = mt2.filter_rows(hl.is_defined(pli[mt2.row_key]))


print("find total n ac5 per cat")
mt2 = mt2.annotate_cols(total_ptv_mac5 = hl.agg.sum(mt2.n_ptv_ac5), 
    total_mpc3_ac5 = hl.agg.sum(mt2.mpc3_ac5), 
    total_mpc2_to_3_ac5 = hl.agg.sum(mt2.mpc2_to_3_ac5), 
    total_mpc_less_2_ac5 = hl.agg.sum(mt2.mpc_less_2_ac5), 
    total_mpc_at_least_2_ac5 = hl.agg.sum(mt2.mpc_at_least_2_ac5), 
    total_syn_ac5 = hl.agg.sum(mt2.syn_ac5))
    # total_n_ac5 = hl.agg.sum(mt2.n_ac5))

print("export")
OUT_FILE_PLI = f'{RESULTS_OUT}/bge_wave1_amr_counts_by_annotation_ac5_in_high_pli_genes.tsv'
mt2.cols().select('total_n_ac5','total_mpc3_ac5','total_mpc2_to_3_ac5','total_mpc_less_2_ac5','total_mpc_at_least_2_ac5','total_syn_ac5','total_ptv_mac5').flatten().export(output=OUT_FILE_PLI)


