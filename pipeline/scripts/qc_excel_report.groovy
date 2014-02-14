// vim: set ts=4:sw=2:expandtab:cindent:
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// VCF to Excel Format Conversion
//
// This script reads the VCF file containing annotations by VEP 
// and produces an Excel file that can be easily viewed and filtered.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

// Quick and simple way to exit with a message
CliBuilder cli = new CliBuilder(usage: "<options> <gatk coverage report prefixes>\n",writer: new PrintWriter(System.err))
err = { msg ->
  System.err.println("\nERROR: " + msg + "\n")
  cli.usage()
  System.err.println()
  System.exit(1)
}

cli.with {
    s "comma separated list of samples to include", longOpt: "samples", args: 1
    t "threshold of coverage for reporting low coverage regions", args:1
    w "threshold of width for reporting low coverage regions (1)", args:1
    o "name of output file", args:1
}


opts = cli.parse(args)
samples = null
if(!opts.s) 
    err "Please provide -s option to specify samples"

println "opts.o = $opts.o"
if(!opts.o) 
    err "Please provide -o option to specify output file name"

int minRegionWidth = 1
if(opts.w)
    minRegionWidth = opts.w.toInteger()

samples = opts.s.split(',')
args = opts.arguments()

// TODO: this is not robust enough sample name matching: 
// need to suffix sample names with underscore
// when the file names are changed to contain library info
matchesSample = { fileName, sample, extension -> new File(fileName).name.startsWith(sample) && fileName.endsWith(extension)}

files = samples.collectEntries { sample ->
        [ 
          sample, 
            [ 
              cov: args.find { matchesSample(it,sample,".sample_cumulative_coverage_proportions") },
              intervals: args.find { matchesSample(it,sample,".sample_interval_statistics") },
              metrics: args.find { matchesSample(it,sample,".metrics") },
              coverage: args.find { matchesSample(it,sample,".cov.txt") }
            ]
        ]
}

// Read the GATK computed coverage levels
covs = samples.collectEntries { sample ->
       if(!files[sample].cov)
           err "Unable to find cumulative coverage file for sample $sample from files in $args"
       lines = new File(files[sample].cov).readLines()*.split('\t')
       [ sample, [ lines[0][1..-1],lines[1][1..-1]*.toFloat() ].transpose().collectEntries()]
}


// Read the Picard deduplication metrics
metrics = samples.collectEntries { sample ->
        if(!files[sample].metrics)
                err "Unable to find Picard metrics file for sample $sample in provided inputs: $args"
        lines = new File(files[sample].metrics).readLines()
        int index = lines.findIndexOf { it.startsWith("LIBRARY") }
        if(index < 0)
            err "Unable to locate LIBRARY line in Picard metrics file ${files[sample].metrics}"

        [ sample, [ [ lines[index].split('\t')[1..-1], lines[index+1].split('\t')[1..-1]*.toFloat() ].transpose().collectEntries() ]]
}

class Block {
    String region
    String chr
    int start
    int end
    String gene
    DescriptiveStatistics stats = new DescriptiveStatistics()

}

// Calculate the contiguous blocks of low coverage from the coverage file
Map sampleBlocks = [:]
Map sampleStats = [:]
Set allGenes = new HashSet()

int threshold = (opts.t == false ? 15 : opts.t.toInteger())
println "Coverage threshold = $threshold"

for(sample in samples) {
    Block block = null
    int lineCount = 0
    int blockCount = 0
    int totalBP = 0
    def coverageStats = new SummaryStatistics()
    def coveragePercentiles = new CoveragePercentile()

    def blocks = []

    def write = {
        blockCount++
        blocks.add(block)
        // println "Writing block ${block.hashCode()} for $gene from $start - $end"
        block = null
    }

    if(!files[sample].coverage)
            err "Unable to find coverage file (*.cov.txt) for sample $sample"

    println "Low cov file for $sample = ${files[sample].coverage}"

    new File(files[sample].coverage).eachLine { line ->
        ++lineCount
        (chr,start,end,gene,offset,cov) = line.split('\t')
        cov = cov.toFloat()
        coverageStats.addValue(cov.toInteger())
        int pos = start.toInteger() + offset.toInteger()
        String region = "$chr:$start"
        ++totalBP
        allGenes.add(gene)

        if(block && block.region != region) 
            write()

        if(cov < threshold) {
            if(!block)  {
               block = new Block(chr:chr, region:region, gene:gene, start:pos)
            }
            block.stats.addValue(cov.toInteger())
            block.end = pos
        }
        else {
            if(block && (block.end - block.start >= minRegionWidth))
                write()
        }

        if(lineCount % 10000 == 0) {
            println(new Date().toString() + "\t" + lineCount + " ($blockCount low coverage blocks observed)")
        }
    }
    sampleBlocks[sample] = blocks
    sampleStats[sample] = [ max: coverageStats.max, min:coverageStats.min, median: coveragePercentiles.getPercentile(50) ]
}


new ExcelBuilder().build {

    // Summary for all samples in the batch
    sheet("Overview") {
        
        // blank row at top
        row {} 

        // sample headings
        bold { row {
                cell("")
                center { cells(samples) }
        }}

        for(metric in ["READ_PAIRS_EXAMINED","UNMAPPED_READS","PERCENT_DUPLICATION"]) {
            row { 
                cell(metric).bold()
                center {
                    cells(samples.collect { metrics[it][metric][0] })
                }
            }
        }

        row {
            cell("Mean Coverage").bold()
            for(s in samples) { cell(sampleStats[s].mean) }
        }

        row {
            cell("Median Coverage").bold()
            center {
                for(s in samples) { cell(sampleStats[s].median) }
            }
        }

        for(depth in [1,10,20,50]) {
            row { center {
                  cell("Frac > ${depth}X").bold()
                  cells(samples.collect { covs[it]["gte_$depth"] })
            }}
        }
    }.autoSize()

    // Per sample summary
    for(sample in samples) {
        println "Writing sheet for $sample"
        sheet(sample) {
            def blocks = sampleBlocks[sample]
            row {
                cell('Total low regions').bold()
                cell(blocks.size())
            }
            row {
                cell('Total low bp').bold()
                cell(blocks.sum { it.end-it.start})
            }
            row {
                cell('Frac low bp').bold()
                if(blocks)
                    cell(blocks.sum { it.end-it.start} / (float)sampleStats[sample].getN())
            }
            row {
                cell('Genes containing low bp').bold()
                cell(blocks*.gene.unique().size())
            }
            row {
                cell('Frac Genes containing low bp').bold()
                cell(blocks*.gene.unique().size() / (float)allGenes.size())
            }

            row {
            }

            bold { row { bottomBorder {
                cells('','','',"Regions","< $threshold x",'','','')
            }}}
            bold { row {
                cells('gene','chr','start','end','min','max','median','length')
            }}

            def lowBed = new File("${sample}.low.bed").newWriter()
            blocks.each { b ->
                b.with {
                    row { 
                        cells(gene, chr, start, end, stats.min, stats.max, stats.getPercentile(50), end-start);
                        cell("ucsc").link("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$chr%3A$start-$end&refGene=pack")
                    }
                    lowBed.println([chr,start,end,stats.getPercentile(50)+'-'+gene].join("\t"))
                }
            }
            lowBed.close()
        }.autoSize()
    }
}.save(opts.o)
