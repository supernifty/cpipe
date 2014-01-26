// vim: set ts=4:sw=4:expandtab:cindent
////////////////////////////////////////////////////////////
//
// Melbourne Genomics Variant Calling Pipeline
//
// This pipeline executes an exome variant calling analysis.
//
// There are a number of software requirements, which you should ensure are 
// satisfied before running the pipeline. These need to be configured 
// in the file called 'config.groovy'. A template is provided in 
// config.groovy.template which you can use to create the file and set 
// the right paths to your tools and reference data.
//
// The best available documentation for setting up the pipeline is located
// here:
//
//   https://sites.google.com/a/student.unimelb.edu.au/melbourne-genomics-pipeline-development/how-tos/setting-up-the-pipeline
// 
// This documentation is temporary and subject to rapid change as the pipeline
// is under heavy development.
//
//    bpipe run -n 4 ../../../pipeline/pipeline.groovy data/*.fastq.gz
// 
// Author: Simon Sadedin, MCRI
//         Members of the Melbourne Genomics
// 
// License: TODO
//
////////////////////////////////////////////////////////////

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
                    "%.gz" * [ fastqc ] + 
                   "L%_R*.gz" * [ align_bwa + index_bam ] +
                   merge_bams +
                   dedup + index_bam + 
                   realignIntervals + realign + index_bam +
                   recal_count + recal + index_bam +
                       [ call_variants, calc_coverage_stats, gatk_depth_of_coverage ]
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

