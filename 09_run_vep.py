# hailctl dataproc start jsealock-schema --vep GRCh38 --requester-pays-allow-all --region us-central1 --project daly-neale-sczmeta --max-idle 15m 
# hailctl dataproc submit jsealock-schema /Users/juliasealock/Desktop/schema/pumas/twist_intervals/09_run_vep.py

import hail as hl

MT= 'gs://schema_jsealock/pumas/twist_intervals/pumas_bge_wave_1_qc_out_scz_and_controls.mt'

ht = hl.read_matrix_table(MT).rows()

# vep GRCh38 has moved to gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json

ht_vep = hl.vep(ht, "gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json")
ht_vep.write("gs://schema_jsealock/pumas/twist_intervals/annotations/vep_annotate_scz_controls_passing_samples.ht", overwrite=True)