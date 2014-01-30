
// vim: ts=4:sw=4:expandtab:cindent:
////////////////////////////////////////////////////////////
//
// Melbourne Genomics Variant Calling Pipeline
//
// Pipeline stage definitions for exome variant calling 
// pipeline. See pipeline.groovy for more information.
//
// Author: Simon Sadedin, MCRI
//         Members of Melbourne Genomics
//
// License: TODO
// 
////////////////////////////////////////////////////////

/**
 * A simple data structure to represent the information
 * associated with a Melbourne Genomics sample.
 */
class SampleInfo {

    /** Sample name */
    String  sample

    /** List of files containing sequence data (FastQ) */
    List    files

    /** Target (flagship) name */
    String  target

    /** List of genes prioritised for the sample */
    List    genes

    /**
     * Parse the given file to extract sample information
     *
     * @return  List of (Map) objects defining properties of each
     *          sample: 
     *              <li>sample name (sample), 
     *              <li>FastQ files (files), 
     *              <li>Name of flagship (target)
     *              <li>Genes to be classed as high priority (genes)
     */
    static parse_sample_info(fileName) {
        def lines = new File(fileName).readLines().grep { !it.trim().startsWith('#') }
        def sample_info = lines.collect { it.split("\t") }.collect { fields ->
                new SampleInfo(
                    sample: fields[0], 
                    files: fields[2].split(",")*.trim().collect {new File(it).parentFile?it:"../data/$it"}, 
                    target: fields[1], 
                    genes:  fields.size()>3?fields[3].split(",")*.trim():[]
                ) 
        }.collectEntries { [it.sample, it] } // Convert to map keyed by sample name
    }
}

def parse_sample_info(fileName) {
    SampleInfo.parse_sample_info(fileName)
}
