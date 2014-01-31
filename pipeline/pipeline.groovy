// vim: set ts=4:sw=4:expandtab:cindent
////////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Variant Calling Pipeline
//
// This pipeline executes the standard variant calling analysis for 
// analysing data for the MGHA project.
//
// The best available documentation for setting up the pipeline is located
// here:
//
//   https://sites.google.com/site/melbournegenomics/how-tos/setting-up-the-pipeline
// 
// Usage:
//
//   bpipe run ../../../pipeline/pipeline.groovy samples.txt
// 
// Author: Simon Sadedin, MCRI
//         Members of the Melbourne Genomics
// 
// Copyright Melbourne Genomics Health Alliance members. All rights reserved.
//
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
//
////////////////////////////////////////////////////////////////////////////

about title: "Melbourne Genomics Demonstration Project Pipeline"

// Create this file by copying config.groovy.template and editing
load 'config.groovy'

// Supporting routines that parse and set up sample information
load 'scripts/Sample.groovy'

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'

inputs "samples.txt" : """
                        File containing two columns, the first being sample name, 
                        the second a comma separated list of FastQ files for the sample
                       """,
       "fastq.gz"    : """
                        Files containing paired end FastQ reads with sequencing data
                       """

sample_info = parse_sample_info(args[0])

// We are specifying that each analysis takes place inside a fixed file structure
// where the parent directory is named according to the batch name. Thus we
// can infer the batch name from the name of the parent directory.
// 
// Note: this variable can be overridden by passing a parameter to bpipe in case
// you are running in a different location.
batch = new File("..").canonicalFile.name

targets = sample_info*.value*.target as Set

run {
    // For each target (flagship) we run the main pipeline in parallel
    targets * [

        set_target_info +
        
        // The first phase is to perform alignment and variant calling for each sample
        sample_info.keySet() * [
               set_sample_info +
                   "%.gz" * [ fastqc ] + check_fastqc +
                   ~"(.*)_R[0-9][_.].*fastq.gz" * [ align_bwa + index_bam ] +
                   merge_bams +
                   dedup + index_bam + 
                   realignIntervals + realign + index_bam +
                   recal_count + recal + index_bam +
                       [ call_variants, calc_coverage_stats, gatk_depth_of_coverage ]
                   + check_coverage
        ] + 

        // The second phase is to merge all the variants for the target (flagship)
        // and then annotate them
        merge_vcf + 
        filter_variants + 
        annotate_vep + index_vcf +
        [ annovar_summarize_refgene + augment_condel + annotate_significance + add_to_database, qc_excel_report]
   ] + 

   // The final phase is to produce the output spreadsheet, 1 per target (flagship)
   targets * [ set_target_info +  vcf_to_excel ]
}

