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
                bold { row {
                    cells('blockid','gene','min','median','chr','start','end','length')
                }}

                // Example row:
                // chr1   247463621   247464578   1   0   1   1
                // chr1    247463621   247464578   2   0   2   1
                // chr1    247463621   247464578   3   0   3   1
                Block block = null
                int threshold =15

                def write = {
                    block.with {
                        row { cells(chr, start, end, stats.min, stats.max, stats.mean, stats.getPercentile(50)) }
                    }
                }
                
                println "Low cov file for $sample = ${files[sample].lowcov}"
                new File(files[sample].lowcov).eachLine { line ->
                    (chr,start,end,gene,offset,cov) = line.split('\t')
                    cov = cov.toFloat()
                    int pos = start.toInteger() + offset.toInteger()
                    String region = "$chr:$start"

                    if(block && block.region != region) {
                        write()
                    }

                    if(cov < threshold) {
                        if(block)
                           block.stats.addValue(cov.toInteger())
                        else
                           block = new Block(region:region, start:pos)
                        block.end = pos
                    }
                    else {
                        if(block) {
                            write()
                            block = null
                        }
                    }
                }
        }
    }
}.save("qc.xlsx")
