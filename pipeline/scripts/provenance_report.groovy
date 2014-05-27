// vim: sw=4:expandtab:ts=4:cindent:
// Start by getting the basic information we will need
si = sample_info[sample]

// Pull out the output files that are for this sample specifically
println "Searching for sample $sample in output files"
sampleFiles = outputGraph.findAllOutputsBy { 
    def branchParts = it.branchPath.split("/")
    return (branchParts.contains(sample)|| branchParts[-1]==si.target) 
}

// Pull out the output files that are for this target region (flagship/cohort)
targetFiles = outputGraph.findAllOutputsBy { it.branchPath.split("/").contains(si.target) }

files = [
    rawbam: sampleFiles.grep { it.stageName == "align_bwa" },
    finalbam: sampleFiles.grep { it.stageName == "recal" },
    vcf: sampleFiles.grep { it.stageName == "call_variants"  && it.outputFile.name.endsWith("vcf") },
    annovarx: sampleFiles.grep { it.stageName == "vcf_to_excel" && it.outputFile.name.endsWith("annovarx.csv") },
    summary: sampleFiles.grep { it.stageName == "summary_pdf" && it.outputFile.name.endsWith(".pdf") }
].collectEntries { key, fs -> [ key, fs.unique { it.outputFile.absolutePath }[0] ] }

tools = sampleFiles.grep { it.tools }                   // Only files with tools
                   .collect { it.tools.split("\n") }    // A single file has multiple tools -> split/ flatten
                   .flatten()*.trim()*.replaceAll(',$','')*.replaceAll('^,','') // remove leading or trailing commas
                   .unique()                            // Don't list tools multiple times
                   .collect { it.split(":") }           // Split the tool:version for each tool 
                   .collectEntries {  [it[0], it[1]] }  // Create a map of tool => version from above split

// Read the version of the pipeline from the version file in the root directory
pipelineVersion = new File(BASE,"version.txt").text.trim()

new PDF().document(outputFile.absolutePath) {

    title "Provenance Report for Study $sample"
    br()

    bold { p "Pipeline" }

    table(cols:2, padding:4) {
        head { cells("Property","Value") }
        cells("Version", pipelineVersion + " / " + new File("revision.txt").text )
        cells("Revision Date", "git log -1 --format=%cd".execute().text)
        cells("Run By", System.properties["user.name"])
        cells("Date", (new Date()).toString())
    }

    br()
    bold { p "Tools" }

    table(cols:2, padding:4, widths: [0.2f,0.4f,0.2f,0.2f]) {
        head {
            cells("Tool","Version")
        }

        tools.each {  tool, version ->
            // BAM files from alignment
            cell tool
            align('center') { cell version }
        }
    }

    br()
    bold { p "Files" }

    table(cols:4,padding:4) {
        head {
            cells("Description","File(s)","Timestamp","File Size")
        }

        fontSize(9) {

            // BAM files from alignment
            cell("Raw Alignment")
            cell(files.rawbam.outputFile.name)
            cell(files.rawbam.timestamp)
            cell(files.rawbam.outputFile.length())

            cell("Final Alignment")
            cell(files.finalbam.outputFile.name)
            cell(files.finalbam.timestamp)
            cell(files.finalbam.outputFile.length())

            cell("Variant Calls")
            cell(files.vcf.outputFile.name)
            cell(files.vcf.timestamp)
            cell(files.vcf.outputFile.length())

            cell("Annotated Variants")
            cell(files.annovarx.outputFile.name)
            cell(files.annovarx.timestamp)
            cell(files.annovarx.outputFile.length())

            cell("Summary PDF")
            cell(files.summary.outputFile.name)
            cell(files.summary.timestamp)
            cell(files.summary.outputFile.length())
        }
    }
}
