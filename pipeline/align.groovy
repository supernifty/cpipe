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
// By default this pipeline will attempt to use all the available cores
// on the computer it runs on. If you don't wish to do that, limit the 
// concurrency by running it with the -n flag:
//
//    bpipe run -n 4 ../../pipeline/pipeline.groovy data/*.fastq.gz
// 
// Assumes: paired end reads
//
// Author: Simon Sadedin, MCRI
//         Members of the Melbourne Genomics Consortium
// 
// License: TODO
//
////////////////////////////////////////////////////////////

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
samples = new File(args[0]).readLines().collectEntries { [ it.split("\t")[0], it.split("\t")[1].split(",") ] } 

run {
    // Align each pair of input files separately in parallel
    samples * [
       set_sample_info +
           "%.gz" * [ fastqc ] +
           alignBWA     
    ] 
}
