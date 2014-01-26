//vim: shiftwidth=4:ts=4:expandtab
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// Script to augment annovar output with a Condel score, until
// Annovar gets around to including it in their annotations.
//
// This script reads the VCF file and Annovar exome summary file
// and then writes out the Annovar file with one extra column, that
// being the Condel score for the Annovar variant.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////

import com.xlson.groovycsv.*
import au.com.bytecode.opencsv.*

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  i 'VCF file to convert to Excel format', args:1
  a 'Annovar file containing annotations', args:1
  o 'Output file (*.csv)', args:1
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
if(!opts.i) 
    err "Please provide -i option to specify VCF file to process"
if(!opts.o) 
    err "Please provide -o option to specify output file name"
if(!opts.a)
    err "Please provide -a option to specify Annovar annotation file"

// First sniff the header line from the Annovar file to get the columns
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

// Not all the fields have headers (why?)
if(!("Qual" in ANNOVAR_FIELDS))
    ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCFIndex vcf = new VCFIndex(opts.i)

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

// These are copied directly from the ##INFO section of an example VCF
// that was processed by VEP. If the flags to VEP are changed, then they
// may need to be updated
VEP_FIELDS = "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|CANONICAL|PolyPhen|SYMBOL|SYMBOL_SOURCE|HGVSc|HGVSp|AA_MAF|EA_MAF|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|PUBMED|ENSP|SIFT|DISTANCE|CLIN_SIG|Condel".split("\\|")

// Output file
def writer = new FileWriter(opts.o)
writer.println(ANNOVAR_FIELDS.join(",")+",Condel")

CSVWriter csvWriter = new CSVWriter(writer);
def annovar_csv = new CsvParser().parse(new File(opts.a).text, separator:',')
int lineIndex = 0
for(av in annovar_csv) {
    ++lineIndex
    if(lineIndex%5000==0)
        println new Date().toString() + "\tProcessed $lineIndex lines"

    def variant = find_vcf_variant(vcf,av,lineIndex)
    if(!variant) {
        println "WARNING: Variant $av.Chr:$av.Start at line $lineIndex could not be found in the original VCF file"
        csvWriter.writeNext((av.values + [""]) as String[])
    }

    // Parse out the VEP annotations - the only one we want is Condel
    def vep = [VEP_FIELDS,variant.info.CSQ.split(",")[0].split("\\|")].transpose().collectEntries()
    
    println "Condel score for $variant.chr:$variant.pos = $vep.Condel"
    csvWriter.writeNext((av.values + [vep.Condel]) as String[])
}
writer.close()
