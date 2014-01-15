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

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'

inputs "samples.txt" : """
                        File containing two columns, the first being sample name, 
                        the second a comma separated list of FastQ files for the sample
                       """,
       "fastq.gz"    : """
                        Files containing paired end FastQ reads with sequencing data
                       """

// Load the samples (todo: replace with something more robust)
// The sample data will come from a real meta data file that is yet to be
// specified.
samples = new File(args[0]).readLines().collectEntries { [ it.split("\t")[0], it.split("\t")[1].split(",") ] } 

// TODO: need to decide how this is specified
flagship = "EPIL"

run {
    // Align the input files for each sample separately in parallel
    samples * [
       set_sample_info +
            "%.gz" * [ fastqc ] +
           align_bwa +
           index_bam +
           dedup + index_bam + 
           realignIntervals + realign + index_bam +
           recal_count + recal + index_bam +
               [ call_variants, calc_coverage_stats, gatk_depth_of_coverage ]
    ] + merge_vcf + annotate + sort_vcf + vcf_to_excel + qc_excel_report
}
