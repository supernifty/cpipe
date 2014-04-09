#!/bin/bash
# vim: ts=4:expandtab:sw=4
########################################################
#
#      Melbourne Genomics Health Alliance
#
########################################################

EMAILS="ssadedin@gmail.com"

#
# Helper functions
function err() {
        (
        echo
        echo "========================= ERROR =================================="
        echo
        echo "$1" | fmt -w 80
        echo
        echo "=================================================================="
        echo
        ) > selftest.error.log

        mail -s "WARNING: Melbourne Genomics SelfTest Failure" $EMAILS  < selftest.error.log

        exit 1
}

function msg() {
        echo
        echo "================================================================"
        echo "$1" | fmt -w 80
        echo "================================================================"
        echo
}

if [ ! -e pipeline ] && [ ! -e designs ];
then
    err "Could not find pipeline directory. Please run this script from the root of the pipeline distribution."
fi

###########################################################
#
# Build a batch directory containing our test data (na18507)
#
if [ -e batches/selftest ];
then
    rm -rf batches/selftest
fi

msg "Creating directories ..."
mkdir -p batches/selftest/data || err "Unable to create selftest batch directory"
mkdir -p batches/selftest/analysis || err "Unable to create selftest batch directory"

msg "Copying test data ..."
cp -v ./batches/na18507/data/000000000_*.fastq.gz batches/selftest/data  || err "Unable to copy selftest data"
cp -v ./batches/na18507/samples.selftest.txt batches/selftest/samples.txt || err "Unable to copy selftest sample file"

pushd batches/selftest/analysis

msg "Running pipeline (should fail due to FastQC) ..."

../../../bpipe run -p VARIANT_DB=variants.selftest.db ../../../pipeline/pipeline.groovy ../samples.txt > output.log

msg "Checking FastQC error detected ..."
[ -e EPIL.xlsx ] && err "Found results spreadsheet but pipeline should have failed due to FastQC failure"

msg "Success: sample 000000000 failed with FastQC error"

msg "Overriding FastQC error ..."
../../../bpipe override check_fastqc.000000000 > check.log 2>&1
grep -q '000000000.*Overridden' check.log || err "Failed to find expected text in check log"

msg "Running pipeline again ..."

../../../bpipe run -p VARIANT_DB=variants.selftest.db ../../../pipeline/pipeline.groovy ../samples.txt > output2.log 2>&1

[ -e EPIL.xlsx ] || err "Failed to find epilepsy result spreadsheet"
[ -e EPIL.qc.xlsx ] || err "Failed to find epilepsy QC  spreadsheet"
[ -e variants/EPIL.merge.filter.vep.vcf ] || err "Failed to find epilepsy variants"

msg "Success: expected output files exist"

EXPECTED_VARIANTS=518

# Count lines in VCF file (this is NOT the number of variants, but a reasonable proxy)
ACTUAL_VARIANTS=`wc variants/EPIL.merge.filter.vep.vcf | awk '{print $1}'`

# TODO: a lot more comparisons to check output, eg: variants
[ "$ACTUAL_VARIANTS" == "$EXPECTED_VARIANTS" ] || \
    err "Epilepsy VCF file has incorrect number of variants (expected=$EXPECTED_VARIANTS, observed=$ACTUAL_VARIANTS)"

msg "Success: correct number of variants observed"

msg "Test succeeded"

mail -s "OK: Melbourne Genomics SelfTest Successful" $EMAILS  <<! Test completed OK 
!

