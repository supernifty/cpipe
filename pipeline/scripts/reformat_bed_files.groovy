
////////////////////////////////////////////////////////////////////////
//
// Conversion script for preparing target region BED files for use in
// the pipeline.
//
// Author:      Simon Sadedin, simon.sadedin@mcri.edu.au
// Date:        February 2014
//
// Copyright Melbourne Genomics Health Alliance members, all rights reserved.
// 
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
//
////////////////////////////////////////////////////////////////////////


BASE="/vlsci/VR0320/shared/production/"

load "../config.groovy"
load "../pipeline_stages_config.groovy"

if(!file(BASE).exists()) {
  println "="*80
  println "Please supply the BASE parameter using -p BASE=... when running the pipeline to point at the pipeline distriubtion root."
  println "="*80
  System.exit(1)
}

flatten = {
    exec """
        echo "Flattening and removing nonstandard chromosomes ..."

        sortBed -i $input.bed | bedtools merge -nms -i - | sed 's/;.*\$//' | grep -v '^chr[0-9]_' > $output.bed
    """
}

annotate = {
    exec """
        echo "Annotating genes ..."

        JAVA_OPTS="-Xmx1g" $BASE/tools/groovy/2.1.1/bin/groovy $BASE/pipeline/scripts/annotate_genes.groovy -r $BASE/tools/annovar/humandb/hg19_refGene.txt $input.bed > $output.bed
    """
}

sort = {
    exec """
        echo "Sorting ..."

        $BASE/tools/IGVTools/2.3.6/igvtools.lowmem sort $input.bed $output.bed
    """
}

rename = {
    produce(branch.name + '.bed') {
        exec "cp $input.bed $output.bed"
    }
}

extract_pgx = {
    msg "Processing PGX variants for target $branch"
    produce(branch.name + '.pgx.vcf') {
        exec """
            cat $DBSNP | JAVA_OPTS="-Xmx1g" groovy 
                -cp $GROOVY_NGS/groovy-ngs-utils.jar
                -e 'pgx = new File("$input.txt").readLines()*.trim(); VCF.filter { it.id in pgx }' > $output.vcf
        """
    }
}

run {
      //"RefSeq_coding_%.bed" * [ flatten + annotate + sort + rename ] +
      "%.pgx.txt" * [ extract_pgx + annotate_vep ] 
    }


