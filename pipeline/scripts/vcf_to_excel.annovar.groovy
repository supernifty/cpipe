// vim: shiftwidth=4:ts=4:expandtab:cindent
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
import com.xlson.groovycsv.*
import au.com.bytecode.opencsv.*
import org.apache.commons.cli.Option

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  s 'comma separated list of samples to include', args:1, required:true
  vcf 'VCF file to convert to Excel format', args: Option.UNLIMITED_VALUES, required:true
  a 'Annovar file containing annotations', args: Option.UNLIMITED_VALUES, required:true
  o 'Name of output file', args:1, required:true
  x 'Comma separated list of functional types to exclude', args:1
  si 'sample meta data file for the pipeline', args:1, required:true
  db 'Sqlite database containing known variants. If known, a column will be populated with the count of times observed.', args:1
  gc 'File listing genes and categories', args:1, required:true
  oocf 'Count at which out-of-cohort variants will cause variants to be filtered from report', args:1, required: true
  pgx 'VCF file containing variants to treat as pharmacogenomic variants (always report)', args:1
  bam 'BAM file for annotating coverage depth where not available from VCF files', args: Option.UNLIMITED_VALUES
  pgxcov 'Coverage threshold below which a pharmocogenomic site is considered untested (15)', args: 1
  annox 'Directory to send Annovar style per-sample summaries to', args: 1, required: true
  log 'Log file for writing information about variants filtered out', args: 1
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

Writer log = new File(opts.log?:'/dev/null').newWriter()

int out_of_cohort_variant_count_threshold = opts.oocf.toInteger() 

def pg_variants = []
if(opts.pgx) 
    pg_variants = VCF.parse(opts.pgx)

int pgx_coverage_threshold = opts.pgxcov ? opts.pgxcov.toInteger() : 15

sample_info = SampleInfo.parse_mg_sample_info(opts.si)
// sample_info = SampleInfo.parse_sample_info(opts.si)

// println "sample_info = $sample_info"

exclude_types = opts.x ? opts.x.split(",") : []

samples = opts.s.split(",")

Map<String,SAM> bams = null
if(opts.bams) {
    bams = opts.bams.collectEntries { def bam = new SAM(it); [ bam.samples[0], bam ] }
    println "="*80
    println "Read ${bams.size()} BAM files for querying read depth:"
    bams.each {  println "Sample: $it.key => $it.value.samFile " }
    println "="*80
}

// Read the header from the first annovar file to find out the column names
ANNOVAR_FIELDS = null
new File(opts.as[0]).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }
// Not all the fields have headers (why?)
if(!("Qual" in ANNOVAR_FIELDS))
    ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Then read all the annovar files
println "Processing ${opts.as.size()} Annovar files"

// connect to database, if specified
sql = null
if(opts.db) {
    Class.forName("org.sqlite.JDBC")
    sql = Sql.newInstance("jdbc:sqlite:${opts.db}")
}

// Read the gene categories
geneCategories = new File(opts.gc).readLines()*.split('\t').collect { [it[0],it[1]] }.collectEntries()

// Order preferred if clinicians need to review output directly
OUTPUT_FIELDS = ["Func", "Gene", "ExonicFunc", "AAChange", "Gene Category", "Priority Index", "Condel", "Conserved", "ESP5400_ALL", "1000g2010nov_ALL", "dbSNP138", "AVSIFT", "LJB_PhyloP", "LJB_PhyloP_Pred", "LJB_SIFT", "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred", "LJB_GERP++", "SegDup", "Chr", "Start", "End", "Ref", "Obs", "Otherinfo", "Qual", "Depth", "#Obs", "RefCount", "AltCount", "CADD"]

OUTPUT_CSV_FIELDS = ["Func","Gene","ExonicFunc","AAChange","Conserved","SegDup","ESP5400_ALL","1000g2010nov_ALL","dbSNP138","AVSIFT","LJB_PhyloP","LJB_PhyloP_Pred","LJB_SIFT","LJB_SIFT_Pred","LJB_PolyPhen2","LJB_PolyPhen2_Pred","LJB_LRT","LJB_LRT_Pred","LJB_MutationTaster","LJB_MutationTaster_Pred","LJB_GERP++","Chr","Start","End","Ref","Obs","Otherinfo","Qual","Depth","Condel","Priority_Index","CADD","Gene Category","Priority Index","CADD","#Obs","RefCount","AltCount"]

CENTERED_COLUMNS = ["Gene Category", "Priority Index", "1000g2010nov_ALL","ESP5400_ALL", "LJB_PhyloP_Pred","LJB_SIFT_Pred","LJB_PolyPhen2","LJB_PolyPhen2_Pred"]

//
// Utility function: Query database to find out how many times a variant has
// been observed a) within the cohort, b) outside the cohort
//
query_variant_counts = { variant, allele, av, studyId ->

    // Look up in database
    def target = sample_info[studyId].target

    // The total observations of the variant
    def variant_count = sql.firstRow(
            """
            select count(*) 
                from variant_observation o, variant v 
                where o.variant_id = v.id and v.chr = $variant.chr and v.pos = $allele.start and v.alt = $allele.alt
            """)[0]

    def variant_count_outside_target = 0
    if(variant_count>1) { // Will always be at least 1 because it will be recorded for the sample we are processing now
        variant_count_outside_target = sql.firstRow("""
                select count(*) 
                from variant_observation o, variant v, sample s
                where o.variant_id = v.id 
                      and v.chr = $variant.chr 
                      and v.pos = $allele.start
                      and v.alt = $allele.alt
                      and o.sample_id = s.id
                      and s.cohort <> $target
                      and s.cohort <> 'AML'
                """)[0] // NOTE: AML excluded explicitly as a special case 
                        // because it includes tumor samples

        log.println "Variant $variant found $variant_count times globally, $variant_count_outside_target outside target $target / AML"
    }

    int variant_count_within_target = sql.firstRow("""
                select count(*) 
                from variant_observation o, variant v, sample s
                where o.variant_id = v.id 
                      and v.chr = $variant.chr 
                      and v.pos = $allele.start
                      and v.alt = $allele.alt
                      and o.sample_id = s.id
                      and s.cohort = $target
                """)[0] 
 
    return [ total: variant_count, other_target: variant_count_outside_target, in_target: variant_count_within_target ]
}

//
// Utility function to collect information about a variant into the columns
// required for export.
//
collectOutputValues = { lineIndex, funcGene, variant, sample, variant_counts, av ->

    // Build up values for the row in a map with the column name as the key
    def outputValues = [:]

    (func,gene) = funcGene
    def aaChange = av.AAChange
    if(gene.indexOf("(")>=0) {
        def geneParts = (gene =~ /(.*)\((.*)\)/)[0]
        gene = geneParts[1].toString()
        aaChange = geneParts[2].toString()
    }
    outputValues.AAChange = aaChange
    outputValues.ExonicFunc = func=="splicing"?"":av.ExonicFunc

    def geneCategory = geneCategories[gene]
    if(sample_info[sample].geneCategories[gene])
        geneCategory = sample_info[sample].geneCategories[gene]

    outputValues["Gene Category"] = geneCategory?:1
    outputValues["Priority Index"] = outputValues["Priority_Index"] = av.Priority_Index
    outputValues["Gene"] = gene
    outputValues["Func"] = func

    for(i in 4..(av.values.size()-3)) {
        outputValues[ANNOVAR_FIELDS[i]] = av.values[i]
    }
    outputValues.CADD = av.columns.CADD != null ? av.CADD: ""

    if(sql) {
        outputValues["#Obs"] = variant_counts.in_target
    }

    outputValues.RefCount=outputValues.AltCount="";
    if(!variant) 
        return outputValues // Cannot annotate allele depths for this variant

    // Try to annotate allele frequencies
    def gt = variant.sampleGenoType(sample)
    if(gt) {
        // Reference depth
        if(gt.containsKey('AD')) {
            outputValues.RefCount=gt.AD[0]

            // Alternate depth depends on which allele
            int altAllele = (variant.alts.size()==1)?1:variant.equalsAnnovar(av.Chr, av.Start.toInteger(), av.Obs)
            outputValues.AltCount = gt.AD[altAllele]
        }
        else {
          System.err.println("WARNING: variant $variant.chr:$variant.pos ($variant.ref/$variant.alt) had no AD info for sample $sample at line $lineIndex")
        }
    }
    else {
      System.err.println("WARNING: variant $variant.chr:$variant.pos ($variant.ref/$variant.alt) had no genotype for sample $sample at line $lineIndex")
    }
    return outputValues
}

//
// Now build our spreadsheet, and export CSV in the same loop
//
try {
    new ExcelBuilder().build {

        for(sample in samples) { // one sample per spreadsheet tab
            def s = sheet(sample) { 
                lineIndex = 0
                sampleCount = 0
                includeCount=0

                // Read the CSV file entirely
                // Sort the annovar output by Priority Index
                String samplePrefix = sample+"."
                String annovarName = opts.as.find{new File(it).name.startsWith(samplePrefix)}
                if(annovarName == null)
                    err "The following samples did not have an Annovar file provided: $sample in Annovar files:\n${opts.as.join('\n')}"

                def annovar_csv = parseCSV(annovarName,',').grep { it.Priority_Index.toInteger()>0 }.sort { -it.Priority_Index.toInteger() }

                // Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
                String vcfName = opts.vcfs.find { new File(it).name.startsWith(samplePrefix) }
                if(vcfName == null)
                    err "The following samples were not found in the VCF file provided: ${sample}"

                VCFIndex vcf = new VCFIndex(vcfName)

                // Write out header row
                bold { row {
                        cells(OUTPUT_FIELDS)
                } }

                println "Priority genes for $sample are ${sample_info[sample].geneCategories.keySet()}"

                // We are going to write out a CSV that is identical to the original annovar output
                // but which includes our custom fields on the end
                // Start by writing the headers
                def writer = new FileWriter("${opts.annox}/${sample}.annovarx.csv")
                writer.println(OUTPUT_CSV_FIELDS.join(","))
                CSVWriter csvWriter = new CSVWriter(writer);
                for(av in annovar_csv) {
                    ++lineIndex
                    if(lineIndex%5000==0)
                        println new Date().toString() + "\tProcessed $lineIndex lines"

                    if(av.ExonicFunc in exclude_types) {
                        log.println "Variant $variant excluded by being an excluded type: $av.ExonicFunc"
                        continue
                    }

                    def variantInfo = vcf.findAnnovarVariant(av.Chr, av.Start, av.End, av.Obs)
                    if(!variantInfo) {
                        println "WARNING: Variant $av.Chr:$av.Start at line $lineIndex could not be found in the original VCF file"
                        log.println "Variant $av.Chr:$av.Start excluded because it could not be identified in the source VCF file"
                        continue
                    }

                    Variant variant = variantInfo.variant
                    if(variant.sampleDosage(sample, variantInfo.allele)==0)
                        continue

                    Map variant_counts = [total: 0, other_target:0]
                    if(sql) {
                        variant_counts = query_variant_counts(variant, variant.alleles[variantInfo.allele], av, sample)
                        if(variant_counts.other_target>out_of_cohort_variant_count_threshold) {
                            log.println "Variant $variant excluded by presence ${variant_counts.other_target} times in other targets"
                            continue
                        }
                    }
                    ++includeCount

                    ++sampleCount
                    # don't split annotations
                    def funcs = [ av.Func ] # av.Func.split(";")
                    def genes = [ av.Gene ] # av.Gene.split(";")

                    [funcs,genes].transpose().each { funcGene ->
                        
                        def outputValues = collectOutputValues(lineIndex, funcGene, variant, sample, variant_counts, av)

                        // Write the row into the spreadsheet
                        row {
                            OUTPUT_FIELDS.each { fieldName ->
                                if(fieldName in CENTERED_COLUMNS) { 
                                    center {cell(outputValues[fieldName])}
                                }
                                else {
                                    cell(outputValues[fieldName]) 
                                }
                            }
                      }

                      // Write Annovar CSV format
                      csvWriter.writeNext( OUTPUT_CSV_FIELDS.collect { fieldName ->
                        outputValues[fieldName] == null ? "" : outputValues[fieldName]
                      } as String[])
                  }
                } // End Annovar variants

                // Now add pharmacogenomic variants
                for(pvx in pg_variants) {

                    values = OUTPUT_FIELDS.collect { "" }

                    // Check if the unfiltered VCF has the variant
                    // Note: there's an issue here about canonicalizing the variant
                    // representation. For now, it's being ignored.
                    def vx = vcf.contains(pvx)

                    List<Map> vepInfos = pvx.vepInfo
                    def genes = vepInfos*.SYMBOL.grep { it != null }.join(",")

                    def state = "Untested"
                    int depth = bams[sample].coverage(pvx.chr, pvx.pos)
                    // println "Queried depth $depth at $pvx.chr:$pvx.pos"
                    // TODO: why is vx sometimes null? should always be genotyped
                    if(vx && depth >= pgx_coverage_threshold) {
                        int allele = vx.findAlleleIndex(pvx.alleles[0])
                        state = vx.sampleDosage(sample, allele) > 0 ? "Present" : "Absent"
                    }
                    else { // Variant not called - but is there coverage?
                        if(depth >= pgx_coverage_threshold)
                            state = "Absent"
                        System.err.println "WARNING : PGX variant $pvx was not genotyped for sample $sample"
                    }

                    // Convert to annovar form since we are using Annovar annotations in the 
                    // rest of the report
                    def annovarVx = vx ? vx.toAnnovar() : pvx.toAnnovar()

                    // Now set the fields that we can
                    def output = [ 
                        'Gene Category': 1, 
                        'Priority Index': 1, 
                        Func: "pharma", 
                        ExonicFunc: state,
                        Gene: genes,
                        dbSNP138: pvx.id,
                        Chr: pvx.chr,
                        Start: pvx.pos,
                        End: pvx.pos + pvx.size(),
                        Ref: annovarVx.ref,
                        Obs: annovarVx.obs,
                        Otherinfo: vx ? (vx.dosages[0] == 1 ? "het" : "hom") : ""
                    ]
                    output.each { k, v ->
                        values[OUTPUT_FIELDS.indexOf(k)] = v 
                    }

                    csv_out = []
                    nvlcell = { cell(it == null ? "" : it ) }
                    row { 
                        OUTPUT_FIELDS.each { fieldName ->
                            //println "Export $fieldName = ${outputValues[fieldName]}"
                            if(fieldName in CENTERED_COLUMNS) { 
                                center { nvlcell(output[fieldName]) } 
                            }
                            else
                            if(fieldName == "ExonicFunc" && state == 'Present') {
                                red {nvlcell(output[fieldName])}  
                            }
                            else {
                                nvlcell(output[fieldName]) 
                            }
                        }
                    }
                    csvWriter.writeNext( ANNOVAR_FIELDS.collect { f -> if(output[f] != null) { String.valueOf(output[f]) } else "" } as String[] )

                }
                csvWriter.close()
            }
            println "Sample $sample has ${sampleCount} / ${includeCount} of included variants"
            try { s/*.autoFilter("A:"+(char)(65+6+samples.size()))*/.autoSize() } catch(e) { println "WARNING: Unable to autosize columns: " + String.valueOf(e) }
            s.setColumnWidth(OUTPUT_FIELDS.indexOf("Gene"),60*256) // 30 chars wide for Gene column
            s.setColumnWidth(OUTPUT_FIELDS.indexOf("AAChange"),30*256) // 60 chars wide for AAChange column
            s.setColumnWidth(OUTPUT_FIELDS.indexOf("Gene Category"),14*256) // 60 chars wide for AAChange column
        }

        sheet("README") {
            row { }
            row { cell("This sheet contains explanations of columns in the previous sheet(s)") }
            row {}
            row { cell("Gene").bold(); cell("The gene affected. A mutation may occur on multiple rows if more than one gene or transcript is affected") }
            row { cell("ESP5400").bold(); cell("Frequency of allele in ESP project (5400 exomes)") }
            row { cell("1000g2010nov_ALL").bold(); cell("Frequency of allele in 1000 Genomes project 2010 Nov release") }
            row { cell("LJB_XXX").bold(); cell("DBNSFP annotations indicating predictions of pathogenicity") }
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
}
finally {
    if(sql)
       sql.close()

    log.close()
}
