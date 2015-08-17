# assume we run from tests

# python unit tests
#python -m unittest discover -s . -p '*_test.py'
python -m unittest discover -s . -p 'update_pipeline_run_id_test.py'

# bpipe tests
../../bpipe run ./dedup_test.groovy ./sample.bam

