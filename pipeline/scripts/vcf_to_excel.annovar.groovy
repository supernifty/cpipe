//vim: shiftwidth=4:ts=4:expandtab
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// VCF to Excel Format Conversion
//
// This script reads the VCF file and Annovar exome summary file
// and produces an Excel file that can be easily viewed and filtered
// by less technical users.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
//
/////////////////////////////////////////////////////////////////////////

import groovy.sql.Sql

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  s 'comma separated list of samples to include', args:1
  i 'VCF file to convert to Excel format', args:1
  a 'Annovar file containing annotations', args:1
  o 'Name of output file', args:1
  x 'Comma separated list of functional types to exclude', args:1
  si 'sample meta data file for the pipeline', args:1
  db 'Sqlite database containing known variants. If known, a column will be populated with the count of times observed.', args:1
}
opts = cli.parse(args)

// Quick and simple way to exit with a message
err = { msg ->
  System.err.println("\nERROR: " + msg + "\n")
  cli.usage()
  System.err.println()
  System.exit(1)
}

if(!opts) {
   cli.usage()
   err "Failed to parse command line options"
}
args = opts.arguments()
if(!opts.s) 
    err "Please provide -s option to specify samples to export"
if(!opts.i) 
    err "Please provide -i option to specify VCF file to process"
if(!opts.o) 
    err "Please provide -o option to specify output file name"
if(!opts.a)
    err "Please provide -a option to specify Annovar annotation file"
if(!opts.si)
    err "Please provide -si option to specify sample meta data file"

sample_info = new Sample().parse_sample_info(opts.si)

println "sample_info = $sample_info"

exclude_types = opts.x ? opts.x.split(",") : []

samples = opts.s.split(",")

// Read all the annovar files
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

// Not all the fields have headers (why?)
if(!("Qual" in ANNOVAR_FIELDS))
    ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCFIndex vcf = new VCFIndex(opts.i)

//println "Loaded ${vcf.variants.size()} variants from VCF file"

missing_samples = samples.grep { !(it in vcf.headerVCF.samples) }
if(missing_samples)
    err "The following samples were not found in the VCF file provided: ${missing_samples.join(',')}"

// connect to database, if specified
sql = null
if(opts.db) {
    Class.forName("org.sqlite.JDBC")
    sql = Sql.newInstance("jdbc:sqlite:${opts.db}")
}

//
// Function to find an Annovar variant in the original VCF file
//
def find_vcf_variant(vcf, av, lineIndex) {
  try {
      int pos = av.Start.toInteger()
      int end = av.End.toInteger()
      return vcf.find(av.Chr,pos-10,end) { variant -> variant.equalsAnnovar(av.Chr,pos,av.Obs) }
  }
  catch(Exception e) {
      try { println "WARNING: unable to locate annovar variant at $lineIndex in VCF ($e)" } catch(Exception e2) { e.printStackTrace() }
  }
}

//
// Now build our spreadsheet
//
new ExcelBuilder().build {

    for(sample in samples) {
        def s = sheet(sample) { 
            lineIndex = 0
            sampleCount = 0
            includeCount=0

            // Read the CSV file entirely
            def annovar_csv = parseCSV(opts.a, ',').grep { it.Priority_Index.toInteger()>0 }.sort { -it.Priority_Index.toInteger() }

            // Write out header row
            bold { row {
                    cells(["Gene Category","Priority Index"] + ANNOVAR_FIELDS[0..-2] + (sql?["#Obs"]:[]) + ['RefCount','AltCount'])
            } }

            println "Priority genes for $sample are ${sample_info[sample].genes}"

            // Sort the annovar output by significance

            for(av in annovar_csv) {
                ++lineIndex
                if(lineIndex%5000==0)
                    println new Date().toString() + "\tProcessed $lineIndex lines"

                if(av.ExonicFunc in exclude_types)
                    continue

                def variant = find_vcf_variant(vcf,av,lineIndex)
                if(!variant) {
                    println "WARNING: Variant $av.Chr:$av.Start at line $lineIndex could not be found in the original VCF file"
                    continue
                }
                ++includeCount
                        
                if(variant.sampleDosage(sample)==0)
                    continue

                ++sampleCount
                def funcs = av.Func.split(";")
                def genes = av.Gene.split(";")
                [funcs,genes].transpose().each { funcGene ->
                    row {
                        (func,gene) = funcGene
                        def aaChange = av.AAChange
                        if(gene.indexOf("(")>=0) {
                            def geneParts = (gene =~ /(.*)\((.*)\)/)[0]
                            gene = geneParts[1].toString()
                            aaChange = geneParts[2].toString()
                        }
                        def exonicFunc = func=="splicing"?"":av.ExonicFunc
                        center {
                            cells(gene in sample_info[sample].genes?1:2, av.Priority_Index)
                        }
                        cells(func,gene,exonicFunc,aaChange)
                        cells(av.values[4..-2])
                    }
                }

              // Look up in database
              if(sql) {
                def variant_count = sql.firstRow("select count(*) from variant_observation o, variant v where o.variant_id = v.id and v.chr = $variant.chr and v.pos = $variant.pos and v.alt = ${av.Obs}")[0]
                cell(variant_count)
              }

              // Try to annotate allele frequencies
              if(variant) {

                  // Reference depth
                  cell(variant.sampleGenoType(sample).AD[0])

                  // Alternate depth depends on which allele
                  int altAllele = (variant.alts.size()==1)?1:variant.equalsAnnovar(av.Chr, av.Start.toInteger(), av.Obs)
                  cell(variant.sampleGenoType(sample).AD[altAllele])
              }
        
              // TODO: check concordance between annoation sources
              /*
              if(av.AA_Match == "False") {
                  excelCells[4].red()
              }
              */
            }
        }
        println "Sample $sample has ${sampleCount} / ${includeCount} of included variants"
        s/*.autoFilter("A:"+(char)(65+6+samples.size()))*/.autoSize()
        s.setColumnWidth(5,60*256) // 30 chars wide for Gene column
        s.setColumnWidth(3,30*256) // 60 chars wide for AAChange column
    }

    sheet("README") {
        row { }
        row { cell("This sheet contains explanations of columns in the previous sheet(s)") }
        row {}
        row { cell("Gene").bold(); cell("The gene affected. A mutation may occur on multiple rows if more than one gene or transcript is affected") }
        row { cell("ESP5400").bold(); cell("Frequency of allele in ESP project (5400 exomes)") }
        row { cell("1000g2010nov_ALL").bold(); cell("Frequency of allele in 1000 Genomes project 2010 Nov release") }
        row { cell("LJB_XXX").bold(); cell("DBNSFP annotations indicating predictions of pathenogenicity") }
        row { cell(""); cell("Numeric: 0 = low impact, 1.0 = high impact") }
        row { cell(""); cell("D=Damaging") }
        row { cell(""); cell("P=Probably Damaging") }
        row { cell(""); cell("T=Tolerated") }
        row { cell(""); cell("N=Neutral") }
        row { cell(""); cell("B=Benign") }
        row { cell(""); cell("See link below for more information") }
        row { cell(""); cell("http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b4.readme.txt").link("http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b4.readme.txt") }
        row {}
        if(opts.x) {
                row{ cell("NOTE:").bold(); cell("The following categories of variant are excluded from this spreadsheet:")}
                row{ cell(""); cell( opts.x ) }
        }
    }.autoSize()
}.save(opts.o)

if(sql)
   sql.close()
