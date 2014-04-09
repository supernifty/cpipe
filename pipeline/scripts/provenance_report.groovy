// vim: sw=4:expandtab:ts=4:cindent:
// Start by getting the basic information we will need
si = sample_info[sample]

// Pull out the output files that are for this sample specifically
sampleFiles = outputGraph.findAllOutputsBy { it.branchPath.split("/").contains { sample } }

// Pull out the output files that are for this target region (flagship/cohort)
targetFiles = outputGraph.findAllOutputsBy { it.branchPath.split("/").contains { si.target } }

files = [
    bam: sampleFiles.grep { it.stageName == "align_bwa" },
    recal: sampleFiles.grep { it.stageName == "recal" },
    realign: sampleFiles.grep { it.stageName == "realign" }
]

new PDF().document(outputFile.absolutePath) {
    title "Provenance Report for Study $sample"

    table {
        head {
            cells("File","Description","Output","Timestamp")
        }


        // BAM files from alignment
        cell "Alignment Files from BWA" 
        cell files.bam*.outputFile*.name.join(",")
        cell files.bam*.timestamp.collect { new Date(it) }*.toString().join(",")

    }
}
