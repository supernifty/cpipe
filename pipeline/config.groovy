// Set location of you reference files here (see below for the files required)
REFBASE="/vlsci/VR0320/shared/hg19"

// Set a good location for storing large temp files here (probably not /tmp)
// This is the scratch directory for Simon's account - not sure if others can access?
TMPDIR="/scratch/VR0193"

// Set location of Picard tools here
PICARD_HOME="/vlsci/VR0320/shared/tools/picard/picard-tools-1.65"

// Set to the reference FASTA file, which must be indexed
REF="$REFBASE/gatk.ucsc.hg19.fasta"

// Set to a VCF file containing DBSNP entries
DBSNP="$REFBASE/dbsnp_132.hg19.vcf"

// Log data from various tools will appear in here
LOG="pipeline.log"

// Set VCF files containing known indels here
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.b37.chr.vcf"

// High confidence SNPs from Illumina Omni chip
SNPS_1000G="$REFBASE/1000G_omni2.5.hg19.sites.vcf"

BASE="/vlsci/VR0320/shared/"

// Location of all the tools that we custom install (eventually
// this may be all of them, for now we are getting some from
// VLSCI)
TOOLS="$BASE/tools"

// Set GATK location here - this pipeline ist tested with 2.3.9
// and will not work with 1.x versions
GATK="$TOOLS/gatk/2.3.9"

// Utilities for making Excel files 
EXCEL="$TOOLS/excel/1.0"

// Location of Bedtools binary
BEDTOOLS="/usr/local/bedtools/2.17.0-intel/bin"

// Utilties for processing NGS data in Java/Groovy
GROOVY_NGS="$TOOLS/groovy-ngs-utils/1.0"

// Set location of Variant Effect Predictor here
// You also need to download the vep_cache data
// and store it in the local directory called 'vep_cache'
// (you can create a symlink to an existing directory with 
// that name if desired).
VEP="$TOOLS/variant_effect_predictor_74"

IGVTOOLS="$TOOLS/IGVTools/2.3.6"

// Location and version of BWA
BWA="/usr/local/bwa/0.7.5a-intel/bin/bwa"

// This is only used for setting read group information in the
// BAM files
PLATFORM="illumina"

// The exome target region to use for the analysis
// Note this is relative to the analysis directory, not
// the directory of this file!
EXOME_TARGET="../design/target_region.bed"

