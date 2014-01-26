//vim: set shiftwidth=4:ts=4:expandtab
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// VCF to Database Script
//
// This script reads the VCF file and Annovar output summary
// and imports the data into a database that can then be used
// to perform structured queries on the data.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           Database driver (Sqlite used by default)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////

import groovy.sql.Sql

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  v 'VCF file to convert to Excel format', args:1
  a 'Annovar file containing annotations', args:1
  b 'Batch name or id', args:1
  db 'Database file to use', args:1
  c 'Create the tables'
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
if(!opts.v) 
    err "Please provide -v option to specify VCF file to process"
if(!opts.db) 
    err "Please provide -db option to specify the database file name"
if(!opts.a)
    err "Please provide -a option to specify Annovar annotation file"
if(!opts.b)
    err "Please provide -b option to specify the batch name"

// Register driver / create database connection
Class.forName("org.sqlite.JDBC")

// Check if we need to create the necessary database tables or not
isNewDb = !new File(opts.db).exists()
sql = Sql.newInstance("jdbc:sqlite:${opts.db}")
if(opts.c || isNewDb) {
  tables = ["""
        create table variant ( 
                id integer primary key asc, 
                chr text,
                pos integer,   
                start integer,
                end integer,
                ref text,
                alt text,
                protein_change text,
                freq_1000g float,
                freq_esp float,
                dbsnp_id text
        );
        CREATE INDEX variant_idx ON variant (chr,pos,alt);
    """,
    """
        create table variant_observation (
               id integer primary key asc, 
               variant_id integer,
               sample_id integer,
               qual float,
               dosage integer,
               date date
        );
        CREATE INDEX variant_observation_idx ON variant_observation (variant_id);
   """,
   """
        create table sample (
               id integer primary key asc, 
               sampleid text,
               batch text
        )
   """
  ]
  tables.each { sql.execute(it) }
  println "Created ${tables.size()} tables"
}

exclude_types = opts.x ? opts.x.split(",") : []

// Read all the annovar fields
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

// Not all the fields have headers (why?)
ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCFIndex vcf = new VCFIndex(opts.v)
samples = vcf.headerVCF.samples

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

for(sample in samples) {

    lineIndex = 0
    sampleCount = 0
    includeCount=0
    annovar_csv = ExcelCategory.parseCSV("", opts.a, ',')
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

        def gt = variant.sampleGenoType(sample)

        ++sampleCount
        def funcs = av.Func.split(";")
        def genes = av.Gene.split(";")
        [funcs,genes].transpose().each { funcGene ->
            (func,gene) = funcGene
            def aaChange = av.AAChange
            if(gene.indexOf("(")>=0) {
                def geneParts = (gene =~ /(.*)\((.*)\)/)[0]
                gene = geneParts[1].toString()
                aaChange = geneParts[2].toString()
            }
            def exonicFunc = func=="splicing"?"":av.ExonicFunc
            // cells(func,gene,exonicFunc,aaChange)
            // cells(av.values[4..-1])
        }


        def sample_row = sql.firstRow("select * from sample where sampleid = $sample")
        if(!sample_row) {
            sql.execute("insert into sample (id, sampleid, batch) values (NULL, $sample, ${opts.b})")
            sample_row = sql.firstRow("select * from sample where sampleid = $sample")
        }

        def variant_row = sql.firstRow("select * from variant where chr=$variant.chr and pos=$variant.pos and alt=$av.Obs")
        if(variant_row) {
                println "Variant $variant.chr:$variant.pos is already known"
        }
        else {
            sql.execute("""
                insert into variant (id,chr,pos,start,end,ref,alt,protein_change,freq_1000g, freq_esp, dbsnp_id) 
                       values (NULL, $variant.chr, $variant.pos, ${av.Start.toInteger()}, ${av.End.toInteger()}, 
                              $av.Ref, $av.Obs, $av.AAChange, 
                              ${av["1000g2010nov_ALL"]},${av["ESP5400_ALL"]}, $av.dbSNP138)
            """)

            variant_row = sql.firstRow("select * from variant where chr=$variant.chr and pos=$variant.pos and alt=$av.Obs")
        }

        sql.execute("""
                insert into variant_observation (id,variant_id,sample_id,qual,dosage, date) 
                            values (NULL, $variant_row.id, $sample_row.id, ${gt.GQ.toDouble()}, ${variant.sampleDosage(sample)}, 'now')
        """)

    }
    println "Sample $sample has ${sampleCount} / ${includeCount} of included variants"
}
sql.close()

