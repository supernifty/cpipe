/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:expandtab:cindent:nowrap */

// Quick and simple way to exit with a message
def err(msg) {
  System.err.println("ERROR: " + msg)
  System.exit(1)
}

// Note: 'AAChange' comes from Annovar but is handled separately
def EXPORTED_ANNOVAR_COLUMNS =
        [ 'Conserved','SegDup','ESP5400_ALL','g1000','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred' ]

def EFFECTS_RANKED = [
    'START_LOST',
    'STOP_GAINED',
    'STOP_LOST',
    'NON_SYNONYMOUS_CODING',
    'EXON',
    'START_GAINED',
    'SYNONYMOUS_CODING',
    'SPLICE_SITE_ACCEPTOR',
    'SPLICE_SITE_DONOR',
    'INTRON',
    '3_PRIME',
    '5_PRIME',
    'DOWNSTREAM',
    'UPSTREAM',
    'INTRAGENIC'
]

def probeColumns(String fileName) {
    // Find the line with the column names and split it
    new File(fileName).withReader { r ->
        String l = r.readLine(),last
        while(l.startsWith('#')) {
            last = l
            l = r.readLine()
        }
        return last.split('\t')
    }
}

def probeInfos(String fileName) {
    // Find the line with the column names and split it
    Map<String,String> infos = [:]
    new File(fileName).withReader { r ->
        String l = r.readLine()
        while(l.startsWith('#')) {
            if(l.startsWith("##INFO"))
                infos[ (l =~ /ID=([A-Z]*)/)[0][1] ] = l

            l = r.readLine()
        }
    }
    return infos
}

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [-s <sample order>] [...csv / vcf files...]")
cli.with {
  s 'define order in which samples should appear in spreadsheet', args:1
  p 'define a list of phenotypically relevant genes to highlight (first column used in tab delimited file)', args:1
  l 'linkify elements (only works for small data sets)'
}
opts = cli.parse(args)
if(!opts) {
   cli.usage()
   err "Failed to parse command line options"
}

args = opts.arguments()

// Genes that will be highlighted in spreadsheet as linked to 
// phenotype of interest
def phenotype_genes = [] as Set
if(opts.p) {
  // PHENOTYPE_GENES = "/home/simons/work/trio/dsdgenes.txt"
  if(!new File(opts.p).exists()) 
        err "Failed to find phenotype gene file $opts.p"
  phenotype_genes = new File(opts.p).readLines().collect { it.split("\t")[0] } as Set
}


if(args.size()<2) {
    System.err.println "Usage: groovy vcf_to_excel.groovy <annovar files> <vcf files>"
    System.exit(1)
}

csvs = args.grep { it.endsWith('.csv') }
vcfs = args.grep { it.endsWith('.vcf') }

if(csvs.size() != vcfs.size()) {
    System.err.println "Please provide the same number of VCF files as CSV files\n"
    System.err.println "Usage: groovy vcf_to_excel.groovy <annovar files> <vcf files>"
    System.exit(1)
}

fileNameParser = new IlluminaFileNameParser()

use(ExcelCategory) {

  def wb = workbook()

  [csvs,vcfs].transpose().each { files ->
    
        annovar_summary = files[0]
        vcf_file = files[1]

        println "Reading Annovar summary $annovar_summary ..."

        // First read the exome annotations from annovar
        def annovar = parseCSV(annovar_summary,
            ['Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','ESP5400_ALL','g1000','dbSNP132','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred','LJB_GERP++','Chr','Start','End','Ref','Obs','Otherinfo'], ',')

        def annovar_variants = [:]
        for(variant in annovar) {
            annovar_variants[variant.Chr + '_' + variant.Start+'_'+variant.Ref+'_'+variant.Obs] = variant
        }

        println "Read ${annovar_variants.size()} annotations"
        def sampleName = "unknown_sample"
        try {
          sampleName = fileNameParser.parse(vcf_file).sample
        }
        catch(Exception e) {
            println "WARNING: Can't decode sample name"
        }

        println "Sample name is $sampleName"

        def sheet = wb.createSheet(sampleName);

        // All lines, including ones we ignored
        int lineCount = 0

        // Count of lines we actually wrote
        int lineIndex = 0

        // Current line being processed in current file
        def currentLine  = "Unknown"

        try {
            // Find the line with the column names. We're doing this to figure out the sample names
            println "Processing $vcf_file"

            def colNames = probeColumns(vcf_file)

            println "$vcf_file has columns $colNames"

            def infoFields = probeInfos(vcf_file)

            println "$vcf_file has info fields ${infoFields.keySet()}"

            def sampleNames = colNames[9..-1]
            def exportSamples = sampleNames
            if(opts.s) {
                exportSamples = opts.s.split(",").toList()
                println "User defined export order: $exportSamples"
                if(!exportSamples.every { it in sampleNames })
                    throw new Exception("Could not find one or more samples $opts.s in sample names from VCF file: $sampleNames")
            }

            println "Exporting samples $exportSamples from $vcf_file"

            // If the file contains CLR (samtools constraint likelihood ratios), add a column for those

            // Write the header line and make it bold
            def exportedColumns = ['gene','chr','start','end','id','rank','effect','qual','depth'] + 
                exportSamples + 
                (exportSamples.size()>1?["scount"]:[]) +
                (exportSamples.size()!=sampleNames.size()?["allcount"]:[]) +
                (infoFields.containsKey("CLR")?["CLR"]:[]) +
                (infoFields.containsKey("SSDNP")?["SSDNP"]:[]) +
                (infoFields.containsKey("FC")?["mut fm cnt"]:[]) +
                (infoFields.containsKey("GC")?["gene fm cnt"]:[]) +
                (phenotype_genes ? ["pheno match"] : []) +
                [ 'Length' ] +
                [ 'DNAChange' ] +
                ['Annotation'] +
                EXPORTED_ANNOVAR_COLUMNS 

            sheet.createRow(lineIndex++).add(exportedColumns) .cellIterator()*.bold()

            VCF.parse(vcf_file) { v ->

              currentLine = v.line // for debugging, save this so the exception can see it

              ++ lineCount

              genes = v.snpEffInfo*.gene
              effs = v.snpEffInfo*.type
              ranks = v.snpEffInfo*.impact

              if(lineCount % 1000 == 0)
                  println lineCount

              String urlPos = v.pos
              String urlChr = URLEncoder.encode(v.chr.replaceAll('chr',''))

              // For trio calling, Samtools outputs CLR in the INFO field.
              // We would like to extract it to a separate column, IF it exists 
              def clr = infoFields.containsKey("CLR") ? v.info.CLR : "0"

              // if(v.info.DP.toInteger()<3)
              //    return

              if(v.qual < 5) 
				return false

              def dosages =exportSamples.collect { v.sampleDosage(it) } 
              if(dosages.every { it == 0 })
                  return false

			  def depths = exportSamples.collect { v.genoTypes[sampleNames.indexOf(it)] }*.DP*.toInteger().collect { it ?: 0 }
 //             println "Depths are $depths"
//System.exit(0)

              if(depths.every { it < 3 })
				return false

			  def gqs = exportSamples.collect { v.genoTypes[sampleNames.indexOf(it)] }*.GQ*.toFloat().collect { it ?: 0.0f }
              if(gqs.every { it < 5 })
				return false

              // There are places where there are overlapping genes
              // To enable easy filtering in the spreadsheet, we actually 
              // write out the same row for each gene
              for(geneAndEffect in [genes,effs,ranks].transpose().unique()) {

                  def (gene,effect,rank) = geneAndEffect
                  //if(!(rank in ['HIGH','MODERATE']))
                  //    continue

                  def row = sheet.row()
                  def genesCell = row.addCell(gene)
                  if(opts.l && v.id == '.') {
                        genesCell.link("http://www.genecards.org/cgi-bin/carddisp.pl?gene=${gene}&search=${gene}")
                        if(phenotype_genes.contains(gene))
                            genesCell.red()
                  }

                  // if(v.id == '.')
                  //    genesCell.link("http://asia.ensembl.org/Homo_sapiens/Gene/Phenotype?g=${URLEncoder.encode(genes[0])};r=$urlChr:$urlPos")

                  // search pubmed
                  // row.addCell(genes.join(',')).link("http://www.ncbi.nlm.nih.gov/pubmed?term=${genes[0]}")

                  row.addCell(v.chr)

                  // UCSC
                  // row.addCell(v.pos).link("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$v.chr%3A$v.pos&dgv=pack&knownGene=pack&omimGene=pack")

                  boolean isDeletion = v.ref.size() > v.alt.size()
                  boolean isInsertion = v.alt.size() > v.ref.size()
                  int len = (v.alt.size() == v.ref.size()) ? 1 : Math.abs(v.alt.size() - v.ref.size())

                  int start = v.pos
                  int end = v.pos
                  if(isDeletion) {
                      start = v.pos+1
                      end = start+len
                  }

                  def posCell = row.addCell(start)

                  // Length
                  row.add(v.pos + len)

                  // ENSEMBL
                  // if(v.id == '.')
                  //    posCell.link("http://asia.ensembl.org/Homo_sapiens/Location/View?g=${URLEncoder.encode(genes[0])};r=$urlChr:$urlPos")

                  row.add(v.id, rank=="HIGH"?2:1, effect, v.qual, v.info.DP)

                  /*
                  row.add(sampleNames.collect { line[it] }.collect { s -> 
                      def alleles = s.split(':')[0].split('/')
                      if(alleles[0] in ['A','C','T','G']) 
                          alleles[0] == alleles[1] ? 2 : 1
                      else
                          alleles.collect { 
                              it == '.' ? -1 : Integer.parseInt(it) 
                          }.countBy { it }[1]?:0
                   })
                   */

                  row.add(dosages)

                  // If there is more than one sample then it's useful to have the count of 
                  // how many total samples the variant was observed in
                  if(exportSamples.size()>1) {
                      row.add(dosages.count { it > 0 })
                  }

                  if(exportSamples.size()!=sampleNames.size()) {
                      row.add(exportSamples.count { v.sampleDosage(it)>0} )
                  }

                  if(infoFields.containsKey("SSDNP")) 
                      row.addCell(v.info.SSDNP)

                  if(infoFields.containsKey("CLR")) {
                      row.addCell(clr?Integer.parseInt(clr):0)
                  }

                  if(infoFields.containsKey("FC")) {
                      row.addCell(v.info.FC)
                  }

                  if(infoFields.containsKey("GC")) {
                      row.addCell(v.info.GC?:"0")
                  }
                  
                  if(phenotype_genes) {
                      row.add(phenotype_genes.contains(gene) ? 1 : 0)
                  }

                  row.add(len) 

                  row.add(v.ref + ' / ' + v.alt)

                  def key = v.chr+'_'+v.pos+'_'+v.ref+'_'+v.alt
                  if(annovar_variants.containsKey(key)) {
                      def variant = annovar_variants[key]
                      row.add(variant.AAChange) //.replaceAll('^.*?:',''))
                      row.add(EXPORTED_ANNOVAR_COLUMNS.collect { variant[it] })
                  }

                  // row.addCell(info.find { it.startsWith('EFF=') }?.substring(4))
                  ++lineIndex
              }
			  return false
            }

            for(int i=0; i<exportedColumns.size(); ++i) {
                sheet.autoSizeColumn((short)i)
            }
            // sheet.setAutoFilter(CellRangeAddress.valueOf("A1:AK"+lineIndex))
        }
        catch(Exception e) {
            //  System.err.println "Failed processing at line $lineCount:\n\n$currentLine\n\n"
            System.err.println "Failed processing at line $lineCount"
            throw e
        }
        println "$lineIndex / $lineCount rows were exported for sample $sampleName"
    }
    def outputFileName = "results.xlsx"
    wb.save(outputFileName)
    println "Output is $outputFileName"
}

