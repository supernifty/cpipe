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

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'

sample_metadata_file = correct_sample_metadata_file( args[0] ) // fix syntax issues and update sample_metadata_file

try {
  sample_info = SampleInfo.parse_mg_sample_info(sample_metadata_file)
}
catch (RuntimeException e) {
  sample_info = SampleInfo.parse_sample_info(sample_metadata_file)
}

// We are specifying that each analysis takes place inside a fixed file structure
// where the parent directory is named according to the batch name. Thus we
// can infer the batch name from the name of the parent directory.
// 
// Note: this variable can be overridden by passing a parameter to bpipe in case
// you are running in a different location.
batch = new File("..").canonicalFile.name

targets = sample_info*.value*.target as Set

samples = sample_info.keySet()

run {
    // Check the basic sample information first
    check_sample_info + check_tools +

    // Create a single BED that contains all the regions we want to call
    // variants in
    create_combined_target + 

    generate_pipeline_id + // make a new pipeline run ID file if required

    // For each target (flagship) we run the main pipeline in parallel
    targets * [

        set_target_info +
        
        // The first phase is to perform alignment and variant calling for each sample
        samples * [
               set_sample_info +
                   "%.gz" * [ fastqc ] + check_fastqc +
                   ~"(.*)_R[0-9][_.].*fastq.gz" * [ align_bwa + index_bam ] +
                   merge_bams +
                   dedup + index_bam + 
                   realignIntervals + realign + index_bam +
                   recal_count + recal + index_bam +
                       [
                         call_variants + call_pgx + merge_pgx +
                            filter_variants + merge_variants +
                            annotate_vep + index_vcf +
                            annovar_summarize_refgene +
                         [add_to_database, augment_condel + annotate_significance, calculate_cadd_scores] + augment_cadd, 
                         calc_coverage_stats + summary_pdf, 
                         gatk_depth_of_coverage 
                       ]
                   + check_coverage
                   + check_karyotype
        ] + qc_excel_report 
   ] + 

   // The 3rd phase is to produce the output spreadsheet, 1 per target (flagship)
   targets * [ set_target_info +  vcf_to_excel ] +

   // Produce a mini bam for each variant to help investigate individual variants
   samples * [ variant_bams ] +

   // And then finally write the provenance report (1 per sample)
   samples * [ provenance_report /* , annovar_to_lovd */ ] +
   
   // And report on similarity between samples
   sample_similarity_report +

   // update metadata and pipeline ID
   create_sample_metadata
}

