# hailctl dataproc start jsealock-schema --project daly-neale-sczmeta --num-workers 20 --num-preemptible-workers 100 --max-idle 15m --packages gnomad
# hailctl dataproc submit jsealock-schema /Users/juliasealock/Desktop/schema/pumas/twist_intervals/10_annotate_vep_on_mt.py
# hailctl dataproc stop jsealock-schema 
import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
from gnomad.utils.vep import process_consequences #see gnmoad python package documentation here: https://broadinstitute.github.io/gnomad_methods/examples/vep.html
hl.plot.output_notebook()

# gsutil -m cp -r \
#   "gs://bd_scz/BGE_Callset_PAISA_QIMR_NeuroMex_KenyaPsych/annotations/fordist_constraint_official_mpc_values_v2_grch38.ht" \
#   "gs://schema_jsealock/pumas/annotations/"

# gsutil -m cp -r \
#   "gs://raw_data_bipolar_dalio_w1_w2/inputs/gnomad.exomes.r2.1.1.non_psych_sites_GRCh38.ht" \
#   "gs://schema_jsealock/pumas/annotations/"


# Input
MT="gs://schema_jsealock/pumas/twist_intervals/pumas_bge_wave_1_qc_out_scz_and_controls.mt"
HT_VEP = "gs://schema_jsealock/pumas/twist_intervals/annotations/vep_annotate_scz_controls_passing_samples.ht"

HT_MPC = 'gs://schema_jsealock/pumas/annotations/fordist_constraint_official_mpc_values_v2_grch38.ht'
HT_GNOMAD_NONPSYCH='gs://schema_jsealock/pumas/annotations/gnomad.exomes.r2.1.1.non_psych_sites_GRCh38.ht'

#HT_CADD = 'gs://hail-datasets-hail-data/CADD.1.4.GRCh38.ht'

# print("load cadd")
# ht_cadd = hl.experimental.load_dataset(name='CADD',
#                                   version='1.4',
#                                   reference_genome='GRCh38',
#                                   region='us',
#                                   cloud='gcp')
#output
HT_ANNOT = 'gs://schema_jsealock/pumas/twist_intervals/annotations/sites_annot.ht'


#####
print('Read in data')
mt = hl.read_matrix_table(MT)
mt.count()
ht_vep = hl.read_table(HT_VEP)
ht_vep.describe()
ht_vep = ht_vep.select('vep')
ht_vep.count()

ht_mpc = hl.read_table(HT_MPC)
#ht_cadd = hl.read_table(HT_CADD)
ht_gnomad_nonpsych = hl.read_table(HT_GNOMAD_NONPSYCH)
# ht_discov = hl.read_table(HT_DISCOVEHR)



print('Annotate with the vep information')
mt = mt.annotate_rows(vep_annot = ht_vep[mt.row_key])
mt = mt.annotate_rows(vep = mt.vep_annot.vep)
mt = mt.drop(mt.vep_annot)
mt.count()

# Need minimal representation of alleles
mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles))
mt = mt.key_rows_by('min_rep')
mt = mt.drop('locus', 'alleles')

mt = mt.annotate_rows(locus = mt.min_rep.locus,
                      alleles = mt.min_rep.alleles)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.drop('min_rep')

print('Count after minimum representation')
mt.count()

# From Duncan: 
# This time, we want to always use the canonical transcript.
# here, use Konrad's function from the collection of gnomad functions to annotate get a 
# one-one correspondance between variant and gene/consequence.
mt = process_consequences(mt)
# this adds worst_consequence_term, worst_csq_for_variant, worst_csq_by_gene and other fields to ds.vep
# and returns an MT with better formatted consequences

# Case builder function by Konrad, modified by Anne:

print("build anno builder")
PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]
MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant",
                 "splice_region_variant"] #Anne: added in splice region variants for now
SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]
OTHER_CSQS = ["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

# Unused: protein_altering_variant, incomplete_terminal_codon_variant (only a tiny number of these; weird variants- Konrad)
# TODO: question, what to do with: "splice_region_variant" (Konrad)

def annotation_case_builder(worst_csq_for_variant_canonical_expr, lof_use_loftee: bool = True, mis_use_polyphen_and_sift: bool = False, mis_use_strict_def: bool = False, syn_use_strict_def: bool = False):
    case = hl.case(missing_false=True)
    if lof_use_loftee:
        case = (case
                .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') #predicted loss-of-function
                .when(worst_csq_for_variant_canonical_expr.lof == 'LC', 'LC'))
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'pLoF')
    if mis_use_polyphen_and_sift:
        case = (case
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence) &
                      (worst_csq_for_variant_canonical_expr.polyphen_prediction == "probably_damaging") &
                      (worst_csq_for_variant_canonical_expr.sift_prediction == "deleterious"), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), "other_missense"))
    else:
        if mis_use_strict_def:
            case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'missense_variant', 'missense')
        else:
            case = case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'missense')
    if syn_use_strict_def:
        case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
    else:
        case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'non_coding')
    return case.or_missing()

mt = mt.annotate_rows(consequence_category = annotation_case_builder(mt.vep.worst_csq_for_variant_canonical, False, True, False, False))


print('Add other annotations')
mt = mt.annotate_rows(mpc = ht_mpc[mt.row_key])
mt = mt.annotate_rows(MPC = mt.mpc.MPC)
mt = mt.drop(mt.mpc)
# mt = mt.annotate_rows(cadd = ht_cadd[mt.row_key])
mt = mt.annotate_rows(inGnomAD_nonpsych = hl.is_defined(ht_gnomad_nonpsych[mt.row_key]))
# mt = mt.annotate_rows(inDiscovEHR = hl.is_defined(ht_discov[mt.row_key].info))

# Annotate whether a site is a SNP, insertion, deletion, or none; and inframe insertion or deletion
mt = mt.annotate_rows(
    type = (hl.case()
    .when((hl.len(mt.alleles[0]) == 1) & (hl.len(mt.alleles[1]) == 1), "SNP")
    .when(hl.len(mt.alleles[0]) < hl.len(mt.alleles[1]), "Insertion")
    .when(hl.len(mt.alleles[0]) > hl.len(mt.alleles[1]), "Deletion")
    .or_missing()),
    infrIndel = ((mt.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_insertion") | (mt.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_deletion")),
    )

# Any variants with len(indel) > 6? 
# pprint(mt.aggregate_rows(hl.agg.counter(mt.alleles.length() > 6))) #None!
mt_rows = mt.rows().repartition(64)
mt_rows.write(HT_ANNOT, overwrite=True)




