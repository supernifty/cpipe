# assume we run from tests

# python unit tests
#python -m unittest discover -s . -p '*_test.py'
python -m unittest discover -s . -p 'update_pipeline_run_id_test.py'
python -m unittest discover -s . -p 'correct_sample_metadata_test.py'
python -m unittest discover -s . -p 'genelist_to_bed_test.py'
python -m unittest discover -s . -p 'find_new_genes_test.py'
python -m unittest discover -s . -p 'update_gene_lists_test.py'

# bpipe tests
#../../bpipe run ./dedup_test.groovy ./sample.bam

