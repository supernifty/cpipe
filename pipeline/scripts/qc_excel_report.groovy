//vim: set ts=4:sw=2:expandtab:cindent
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


CliBuilder cli = new CliBuilder(usage: "<options> <gatk coverage report prefixes>")
cli.with {
    s "comma separated list of samples to include", longOpt: "samples", args: 1
}

opts = cli.parse(args)

samples = null
if(opts.s) {
    samples = opts.s.split(',')
}

args = opts.arguments()

err = { msg ->
  System.err.println "\nERROR: $msg\n"
  System.exit(1)
}

// The files are named similarly to this:
// NA18507.dedup.realign.recal.cov.sample_interval_statistics

//println "args = $args"

// TODO: need to suffix sample names with underscore
// when the file names are changed to contain library info
files = samples.collectEntries { sample ->
        [ 
          sample, 
            [ 
              cov: args.find { it.startsWith(sample) && it.endsWith(".sample_cumulative_coverage_proportions")},
              intervals: args.find { it.startsWith(sample) && it.endsWith(".sample_interval_statistics")},
              metrics: args.find { it.startsWith(sample) && it.endsWith(".metrics") },
              lowcov: args.find { new File(it).name.startsWith(sample) && it.endsWith(".cov.txt") }
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
        lines = new File(files[sample].metrics).readLines()
        int index = lines.findIndexOf { it.startsWith("LIBRARY") }
        if(index < 0)
            err "Unable to locate LIBRARY line in Picard metrics file"

        [ sample, [ [ lines[index].split('\t'), lines[index+1].split('\t')*.toFloat() ].transpose().collectEntries() ]]
}

class Block {
    String region
    String chr
    int start
    int end
    String gene
    DescriptiveStatistics stats = new DescriptiveStatistics()

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
        

        for(depth in [1,10,20,50]) {
            row { center {
                  cell("Frac > $depth").bold()
                  cells(samples.collect { covs[it]["gte_$depth"] })
            }}
        }
    }.autoSize()

    // Per sample summary
    for(sample in samples) {
        sheet(sample) {

                Block block = null
                int threshold =15
                int lineCount = 0
                int blockCount = 0
                int totalBP = 0
                Set allGenes = new HashSet()

                def blocks = []

                def write = {
                    blockCount++
                    blocks.add(block)
                    // println "Writing block ${block.hashCode()} for $gene from $start - $end"
                    block = null
                }
                
                println "Low cov file for $sample = ${files[sample].lowcov}"

                new File(files[sample].lowcov).eachLine { line ->
                    ++lineCount
                    (chr,start,end,gene,offset,cov) = line.split('\t')
                    cov = cov.toFloat()
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
                        if(block) 
                            write()
                    }

                    if(lineCount % 10000 == 0) {
                        println(new Date().toString() + "\t" + lineCount + " ($blockCount low coverage blocks observed)")
                    }
                }

                row {
                    cell('Total low regions').bold()
                    cell(blocks.size())
                }
                row {
                    cell('Total low bp').bold()
                    cell(blocks.sum { it.end-it.start})
                }
                row {
                    cell('Percent low bp').bold()
                    cell(blocks.sum { it.end-it.start} / (float)totalBP)
                }
                row {
                    cell('Genes containing low bp').bold()
                    cell(blocks*.gene.unique().size())
                }
                row {
                    cell('Percent Genes containing low bp').bold()
                    cell(blocks*.gene.unique().size() / (float)allGenes.size())
                }

                row {
                }

                bold { row {
                    cells('gene','chr','start','end','min','max','median','length')
                }}

                def lowBed = new File("${sample}.low.bed").newWriter()
                blocks.each { b ->
                    b.with {
                        row { cells(gene, chr, start, end, stats.min, stats.max, stats.getPercentile(50), end-start) }
                        lowBed.println([chr,start,end,stats.getPercentile(50)+'-'+gene].join("\t"))
                    }
                }
                lowBed.close()
        }.autoSize()
    }
}.save("qc.xlsx")
