#!/bin/bash
#vim: ts=4:expandtab:sw=4:
############################################################
#
# Installation check script for Melbourne Genomics Pipeline
# dependencies. 
#
# This script does not do the complete job, but it checks
# a number of things and can improve over time.
#
# Author: Simon Sadedin, ssadedin@mcri.edu.au
#
############################################################

# Default versions of tools. Would be nice to share this with
# config.groovy
: ${BWA="0.7.5a"}
: ${GATK="2.8.1-g932cd3a"}
: ${SAMTOOLS="0.1.19"}
: ${BEDTOOLS="2.18.2"}
: ${ANNOVAR:="2013aug23"}
: ${VEP="74"}

# Helper functions
function err() {
        echo
        echo "========================= ERROR =================================="
        echo
        echo "$1" | fmt
        echo
        echo "=================================================================="
        echo
        exit 1
}

function warn() {
        echo
        echo "============================================================"
        echo "WARNING: $1" | fmt
        echo "============================================================"
        echo
}


function msg() {
        echo
        echo "============================================================"
        echo "$1"
        echo "============================================================"
        echo
}

function compile() {
    PROGRAM="$1"
    msg "Check $PROGRAM is compiled ..."
    if [ ! -e $PROGRAM ];
    then
            read -p "$PROGRAM does not seem to be compiled: do you want me to compile it? (y/n): "
            if [ "$REPLY" == "y" ];
            then
                pushd `dirname $PROGRAM`
                [ ! -e Makefile ] && cd ..
                [ ! -e Makefile ] && err "Cannot find Makefile for $PROGRAM"
                make || err "$PROGRAM failed to compile"
                popd
                [ ! -e $PROGRAM ] && err "Could not find $PROGRAM even after compiling it"
            else
                warn "You will need to compile $PROGRAM or edit config.groovy to point to your installation manually"
            fi
    fi
}

eval `sed 's/\/\/.*$//' pipeline/config.groovy` 

[ ! -e pipeline/config.groovy ] && \
        err "You haven't created the file pipeline/config.groovy yet: you need to copy the file pipeline/config.groovy.template and edit it."

msg "Checking dependencies ..."

[ -e pipeline ] || \
        err "I can't see the pipeline directory. Maybe you didn't clone the repository, or you're running this script from the wrong location?"

compile "$BWA"

compile "$SAMTOOLS/samtools"

compile "$BEDTOOLS/bin/bedtools"

msg "Check GATK is downloaded and available"
[ -e $GATK/GenomeAnalysisTK.jar ] || \
        err "Could not locate GATK jar file. Please download and install GATK to tools/gatk/$GATK/" 

msg "Check Annovar is downloaded and available"
[ -e $ANNOVAR/annotate_variation.pl ] || \
        err "Could not locate Annovar script. Please download and install Annovar to tools/annovar/$ANNOVAR/"

msg "Check Annovar database exists"
for i in hg19_snp138.txt hg19_avsift.txt hg19_esp5400_all.txt hg19_refGene.txt hg19_ALL.sites.2010_11.txt hg19_phastConsElements46way.txt; 
do
    [ -e tools/annovar/humandb/$i ] || {
        err "Failed to find all necessary Annovar data files ($i): please use Annovar downdb to download all the data files. See pipeline/scripts/pipeline/scripts/download_annovar_db.sh for assistance."
    }
done

msg "Check VEP database downloaded for version $VEP_VERSION..."
[ -e tools/vep/vep_cache/homo_sapiens/$VEP_VERSION/1 ] || \
    err "Failed to find downloaded VEP data files. Please install homosapiens cache and FASTA files using: 'cd $VEP; perl INSTALL.pl -c ../vep_cache"

msg "Check that reference FASTA exists"
[ -e "$REF" ] || err "Reference FASTA file $REF could not be found. Please place it there or change config.groovy to point to your reference"

[ -e "$REF.fai" ] || err "Reference FASTA file $REF is not indexed. Please run samtools faidx to index it"

[ -e "$DBSNP" ] || err "The DBSNP file $DBSNP does not exist. Please download it."

[ -e "$GOLD_STANDARD_INDELS" ] || err "The indel file $GOLD_STANDARD_INDELS does not exist. Please download it from the GATK resource bundle."

msg "Success: all the dependencies I know how to check are OK"

