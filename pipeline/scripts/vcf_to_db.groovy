//vim: shiftwidth=4:ts=4:expandtab:
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
  v 'VCF file to convert to Excel format', args:1, required: true
  a 'Annovar file containing annotations', args:1, required: true
  b 'Batch name or id', args:1, required: true
  cohort 'Cohort / target /flagship name', args:1, required: true
  db 'Database file to use', args:1, required: true
  idmask 'Regular expression used to mask sample ids to allow multiple ids to match to single individuals', args:1
  c 'Create the tables'
}
opts = cli.parse(args)

// Quick and simple way to exit with a message
err = { msg ->
  System.err.println("\nERROR: " + msg + "\n")
  System.err.println()
  System.exit(1)
}

if(!opts) {
   err "Failed to parse command line options"
}

args = opts.arguments()

// Register driver / create database connection
Class.forName("org.sqlite.JDBC")

// Check if we need to create the necessary database tables or not
isNewDb = !new File(opts.db).exists()
sql = Sql.newInstance("jdbc:sqlite:${opts.db}")
if(opts.c || isNewDb) {
  tables = ["""
        create table variant ( 
                id integer primary key asc, 
                chr text NOT NULL,
                pos integer NOT NULL,   -- the position as specified in VCF
                start integer,          -- the start position as specified by Annovar (different to VCF)
                end integer,            -- end position as specified by Annovar
                ref text NOT NULL,      -- reference sequence
                alt text NOT NULL,      -- alternate (observed) sequence
                protein_change text,    -- protein change (where applicable)
                freq_1000g float,       -- frequency in 1000 genomes (if known)
                freq_esp float,         -- frequency in ESP (if known)
                dbsnp_id text           -- DBSNP ID (if known)
        );
    """,
    """
        CREATE INDEX variant_idx ON variant (chr,pos,alt);
    """,
    """
        create table variant_observation (
               id integer primary key asc, 
               variant_id integer references variant(id),
               sample_id integer references sample(id),
               qual float,
               dosage integer,  -- how many copies of the variant (1=het, 2=hom)
               created datetime NOT NULL,
               UNIQUE (variant_id, sample_id) ON CONFLICT ROLLBACK
        );
   """,
   """
        CREATE INDEX variant_observation_idx ON variant_observation (variant_id);
   """,
   """
        create table sample (
               id integer primary key asc, 
               study_id text NOT NULL,
               batch text,
               cohort text,
               created datetime NOT NULL,
               UNIQUE (study_id) ON CONFLICT ROLLBACK
        );
   """,
   """
        CREATE INDEX sample_index ON sample (study_id);
   """
  ]
  tables.each { sql.execute(it) }
  println "Created ${tables.size()} tables / indexes "
}

if(sql) {
    sql.execute("PRAGMA synchronous = 0;")
    sql.execute("PRAGMA journal_mode = WAL;")
    sql.execute("BEGIN TRANSACTION;")
}

// Read all the annovar fields
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

// Not all the fields have headers (why?)
ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCFIndex vcf = new VCFIndex(opts.v)
samples = vcf.headerVCF.samples

add_variant_to_db = { variant, allele, av, dosage, sample_row ->

    def variant_row = sql.firstRow("select * from variant where chr=$variant.chr and pos=$allele.start and alt=$allele.alt")
    if(variant_row == null) {
        sql.execute(""" insert into variant (id,chr,pos,start,end,ref,alt,protein_change,freq_1000g, freq_esp, dbsnp_id) values (NULL, $variant.chr, $allele.start, ${av?.Start?.toInteger()}, ${av?.End?.toInteger()}, ${variant.ref}, $allele.alt, ${av?.AAChange}, ${av?av["1000g2010nov_ALL"]:null},${av?av["ESP5400_ALL"]:null}, ${av?.dbSNP138}) """)
        variant_row = sql.firstRow("select * from variant where chr=$variant.chr and pos=$allele.start and alt=$allele.alt")
    }

    def gt = variant.sampleGenoType(sample_row.study_id)

    def variant_obs = sql.firstRow("select * from variant_observation where sample_id = ${sample_row.id} and variant_id = ${variant_row.id}")
    if(!variant_obs) {
        sql.execute("""
                insert into variant_observation (id,variant_id,sample_id,qual,dosage, created) 
                            values (NULL, $variant_row.id, ${sample_row.id}, ${gt?.GQ?.toDouble()}, ${dosage}, datetime('now'))
        """)
        return false // variant did not already exist
    }
    else
        return true // variant already existed in database
}

sample_id_mask = opts.idmask ? opts.idmask : false

for(studyId in samples) {

    // The study ID is the unmasked form, which is what is in the VCF file
    // For entry into the database we will apply the mask to the study id
    def sample = sample_id_mask ? (studyId =~ sample_id_mask)[0] : studyId

    def sample_row = sql.firstRow("select * from sample where study_id = $sample")
    if(!sample_row) {
        sample_row = sql.execute("""
            insert into sample (id, study_id, batch, cohort, created) values (NULL, $sample, ${opts.b}, ${opts.cohort}, datetime('now'))
        """)
        sample_row = sql.firstRow("select * from sample where study_id = $sample")
    }

    lineIndex = 0
    sampleCount = 0
    includeCount=0
    existingCount=0
    annovar_csv = ExcelCategory.parseCSV("", opts.a, ',')
    ProgressCounter.withProgress { 
        for(av in annovar_csv) {
            ++lineIndex
            if(lineIndex%5000==0)
                println new Date().toString() + "\tProcessed $lineIndex lines"

            def variantInfo = vcf.findAnnovarVariant(av.Chr, av.Start, av.End, av.Obs)
            if(!variantInfo) {
                println "WARNING: Variant $av.Chr:$av.Start at line $lineIndex could not be found in the original VCF file"
                continue
            }

            ++includeCount

            Variant variant = variantInfo.variant
            int dosage = variant.sampleDosage(studyId, variantInfo.allele)
            if(dosage==0)  // Sample doesn't have the variant: some other sample in the cohort must have it
                continue

           ++sampleCount
            if(add_variant_to_db(variant,variant.alleles[variantInfo.allele],av,dosage, sample_row)) {
                existingCount++
            }
            count()
        }
    }

    println "=" * 80
    println "Now including variants not annotated by Annovar"
    
    VCF vcfFile = VCF.parse(opts.v) { variant ->
        variant.alleles.eachWithIndex { allele, index ->
            int dosage = variant.sampleDosage(studyId, index)
            add_variant_to_db(variant, allele, null, dosage, sample_row)
        }
    }

    println "Sample $studyId has ${sampleCount} / ${includeCount} of included Annovar variants"
    if(existingCount>0) {
        println "WARNING: sample $studyId has $existingCount variants that were already registered in the database"
        println "WARNING: this sample may have been re-processed or re-sequenced."
    }
}
sql.execute("COMMIT;")
sql.close()

