// vim: shiftwidth=4:ts=4:expandtab:
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// Sequencing Summary Report
//
// This script reads the coverage and sample information and
// creates a summary report suitable for high level evaluation
// of the sequencing quality and success in PDF form.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////
import java.awt.Color

/////////////////////////////////////////////////////////////////////////
//
// Parse command line options
//
/////////////////////////////////////////////////////////////////////////
CliBuilder cli = new CliBuilder(usage: "qc_pdf.groovy <options>")

cli.with { 
    cov "Coverage file from Bedtools", args:1, required:true
    study "ID of study for which this QC report is being generated", args:1, required: true
    meta "Meta data for the sample to produce the QC report for", args:1, required:true
    threshold "Coverage threshold for signalling a region as having satisfactory coverage", args:1, required: true
    classes "Percentages of bases and corresponding classes of regions in the form: GOOD:95:GREEN,PASS:80:ORANGE,FAIL:0:RED", args:1, required:true
    o    "Output file name (PDF format)", args: 1, required:true
}

opts = cli.parse(args)
if(!opts)
  System.exit(1)

err = { System.err.println "\nERROR: $it\n"; System.exit(1) }
if(!new File(opts.cov).exists()) 
  err "The input coverage file ($opts.cov) does not exist or could not be accessed"

if(!opts.threshold.isInteger())
  err "Please provide the coverage threshold as an integer between 1 and 1000"

int coverageThreshold = opts.threshold.toInteger()

if(!new File(opts.meta).exists())
  err "The sample meta data file ($opts.meta) does not exist or could not be accessed"

Map samples = SampleInfo.parse_sample_info(opts.meta)
if(!samples.containsKey(opts.study))
  err "The provided meta data file ($opts.meta) did not contain meta information for study $opts.study"

SampleInfo meta = samples[opts.study]

/////////////////////////////////////////////////////////////////////////
//
// Parse user specified classes into colors and levels
//
/////////////////////////////////////////////////////////////////////////
def classes
try {
  classes = opts.classes.split(",")*.trim().collect { it.split(":") }.collect { 
    [it[0], it[1].toInteger(),it[2]] 
  }.sort{ -it[1] }
}
catch(Exception e) {
  err "The class string $opts.classes couldn't be parsed. Please use the format: <class1>:<percentage>:<color>,<class2>:<percentage>..."
}

println "Percentage thresholds are: $classes"

println "Meta info = $meta"

/////////////////////////////////////////////////////////////////////////
//
// Calculate Statistics
//
/////////////////////////////////////////////////////////////////////////    
CoverageStats stats = null
int totalBP = 0
int totalOK = 0
String currentGene = null

List<Map> geneReport = []

CoverageStats allGeneStats = new CoverageStats(1000)

ProgressCounter.withProgress {
  new File(opts.cov).eachLine { line ->
    def (chr,start,end,gene,offset,cov) = line.split("\t")
    (start,end,offset,cov) = [start,end,offset,cov]*.toInteger()
    if(gene != currentGene) {
      if(stats != null) {
        def geneSummary = [
              gene: gene, 
              fracOK: totalOK / (float)totalBP,
              median: stats.median
            ]

        geneReport.add(geneSummary)
      }
      stats = new CoverageStats(1000) 
      currentGene = gene
    }
    if(cov > coverageThreshold) 
      ++totalOK

    stats.addValue(cov)
    allGeneStats.addValue(cov)

    ++totalBP
    count()
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Generate PDF
//
/////////////////////////////////////////////////////////////////////////
// Read the coverage file and write out a line in the PDF for each gene containing how many 
// bp are below the threshold
new PDF().document(opts.o) {

  title("Sequencing Summary Report for Study $opts.study")


  bold { p("Summary Data") }

  def hd = { text -> bg("#eeeeee") { bold { cell(text) } } }
  table(cols:2,padding:4) {
    hd("Batch"); cell(meta.batch);
    hd("Study ID"); cell(meta.sample);
    hd("Sex"); cell(meta.sex);
    hd("Disease Cohort"); cell(meta.target);
    hd("Hospital / Institution"); cell(meta.institution);
    hd("Ethnicity"); cell(meta.ethnicity);
    hd("Prioritzed Genes"); cell(meta.geneCategories.collect { it.key }.join(","));
    hd("Consanguinity Status"); cell(meta.consanguinity.name().replaceAll("_"," "));
    hd("Sample Type (tumor/normal)"); cell(meta.sampleType);
    hd("Sequencing Dates"); fontSize(10) { cell(meta.sequencingDates*.format("yyyy-MM-dd")?.unique()?.join(", ")); }
    hd("DNA Collection Dates"); fontSize(10) { cell(meta.dnaDates*.format("yyyy-MM-dd")?.unique()?.join(", ")); }
    hd("Sequencing Machines"); cell(meta.machineIds?.unique()?.join(","));
  }
  br()

  bold { p("Coverage Summary") }
  table(cols:2,padding:4) {
    hd("Reported Mean Coverage"); cell(meta.meanCoverage);
    hd("Observed Mean Coverage"); cell(String.format("%2.1f",allGeneStats.mean));
    hd("Observed Median Coverage"); cell(allGeneStats.median);
  }

  br()

  bold { p("Gene Summary") }
  table {
    bg("#eeeeee") { head {
      cells("Gene", "Perc > ${coverageThreshold}X","Median", "OK?")
    } }

    for(geneSummary in geneReport) {
        cell(geneSummary.gene)
        align("center") {
          cells(String.format("%2.1f%%",100*geneSummary.fracOK), geneSummary.median)
        }

        // Color depends on class
        def clazz = classes.find { geneSummary.fracOK*100 > it[1] }
        color(clazz[2]) {
          cell(clazz[0])
        }
    }
  }
}
