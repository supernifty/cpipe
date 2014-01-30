// vim: cindent:ts=4:expandtab:sw=2
/**
 * Add genes to a bed file that are annotated in a second bed file,
 * or in RefSeq table downloaded from UCSC. 
 */
CliBuilder cli = new CliBuilder(usage:"annotate_genes.groovy [-r <ucsc refgene annotations>] [-c] <bed_to_annotate> <bed_with_genes>")

cli.with {
    r "UCSC RefGene annotation file", args:1
    c "Check mode: do not write output BED, only check inputs and return status 1 if errors found"
}
opts = cli.parse(args)
args = opts.arguments()

err = System.err
if(args.size()<1) {
        err.println("Usage: annotate_genes.groovy -r <annotations> <bed_to_annotate> <bed_with_genes>\n")
        err.println("Adds genes to column 4 of a BED file based on gene annotations from a second bed file")
        System.exit(1)
}

geneBed = null
if(args.size() > 1) {
    geneBed = new BED(withExtra:true, args[1])
    geneBed.withExtra = true
    geneBed.load()
}

// Load refgene database to check annotations here
refGenes = null
if(opts.r) {
    err.println "Loading genes from RefGene database ..."
    refGenes = new BED()
    new File(opts.r).eachLine { line->
        def fields = line.split('\t')
        def (chr,start,end,gene) = [fields[2],fields[4].toInteger(),fields[5].toInteger(),fields[12]]
        
        refGenes.add(chr,start,end,gene)
    }
    err.println "Finished."
}

// Let's also check that each gene is correctly annotated to only
// a single chromosome, since Agilent seems to have huge problems
// with that
gene_chr = [:]

unknownCount = 1
misannotatedSize = 0
introducedGenes = [] as Set

new BED(args[0]).eachRange { chr, start,end ->
        // Find overlapping intervals in gene bed
        def overlaps = []
        if(overlaps) {
            overlaps = geneBed.getOverlaps(chr,start,end)
        }
        def gene = overlaps ? overlaps[0].extra : "Unknown"
        if(!gene_chr[gene])
                gene_chr[gene] = chr

        if(refGenes != null) {
            overlappingGenes = refGenes.getOverlaps(chr,start,end)*.extra.unique()
            if(!(gene in overlappingGenes)) {
                if(geneBed)
                    err.println "ERROR: Gene $gene is annotated to both $chr and ${gene_chr[gene]} at $chr:$start-${end}. Correct gene is ${overlappingGenes} from UCSC annotations"
                misannotatedSize += (end-start)
                // In case there are multple overlapping genes, we choose the shortest one.
                // the reason is the cases where there is an overlap is usually some weird thing
                // with a long name (eg: LOC100505826) and the short name is what people actually want
                gene = overlappingGenes.min { it.size() }
                if(gene == null) {
                    gene = "Intergenic_${unknownCount++}"
                }
                else {
                    introducedGenes << gene
                }
            }
        }
        else {
            if(gene_chr[gene] != chr) {
                err.println "ERROR: Gene $gene is incorrectly annotated to both chromosome $chr and ${gene_chr[gene]}"
            }
        }
        // err.println "(aborted)"
        // System.exit(1)
        if(!opts.c)
            println "$chr\t$start\t$end\t$gene"
}

if(misannotatedSize>0) {
    err.println "WARNING: regions annotated incorrectly - total size: $misannotatedSize"
}

if(introducedGenes.size()>0) {
    err.println "WARNING: The following ${introducedGenes.size()} genes were added based on UCSC annotations:\n$introducedGenes"
}

if(unknownCount>1) {
    err.println "WARNING: ${unknownCount-1} intergenic regions were found" 
}

if(opts.c && (misannotatedSize || introducedGenes || unknownCount)) {
        System.exit(1);
}
