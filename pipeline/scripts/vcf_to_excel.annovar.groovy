//vim: set shiftwidth=4:ts=4:expandtab
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// VCF to Excel Format Conversion
//
// This script reads the VCF file containing annotations by VEP 
// and produces an Excel file that can be easily viewed and filtered
// by less technical users.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////


// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  s 'comma separated list of samples to include', args:1
  i 'VCF file to convert to Excel format', args:1
  a 'Annovar file containing annotations', args:1
  t 'Name for main spreadsheet tab (should reflect batch, sample group, etc)', args:1
  o 'Name of output file', args:1
  x 'Comma separated list of functional types to exclude', args:1
}
opts = cli.parse(args)
if(!opts) {
   cli.usage()
   err "Failed to parse command line options"
}

// Quick and simple way to exit with a message
err = { msg ->
  System.err.println("\nERROR: " + msg + "\n")
  cli.usage()
  System.err.println()
  System.exit(1)
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
if(!opts.t)
    err "Please provide -t option to specify title for spreadsheet"

exclude_types = opts.x ? opts.x.split(",") : []

samples = opts.s.split(",")

// Read all the annovar files
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

// Not all the fields have headers (why?)
ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCFIndex vcf = new VCFIndex(opts.i)

//println "Loaded ${vcf.variants.size()} variants from VCF file"

missing_samples = samples.grep { !(it in vcf.headerVCF.samples) }
if(missing_samples)
    err "The following samples were not found in the VCF file provided: ${missing_samples.join(',')}"

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
            annovar_csv = parseCSV(opts.a, ',')
            for(av in annovar_csv) {
                if(lineIndex == 0) {
                    // Write out the header columns
                    bold {
                        row {
                            cells(ANNOVAR_FIELDS + ['RefCount','AltCount'])
                            cells(samples)
                        }
                    }
                }
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
                        cells(func,gene,exonicFunc,aaChange)
                        cells(av.values[4..-1])
                    }
                }

              // Try to annotate allele frequencies
              if(variant)
                  cells(variant.sampleGenoType(sample).AD)

        
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
        s.setColumnWidth(3,30*256) // 30 chars wide for Gene column
        s.setColumnWidth(1,60*256) // 60 chars wide for AAChange column
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
