// vim: ts=4:sw=4:expandtab:cindent:number
////////////////////////////////////////////////////////////
//
// Melbourne Genomics Variant Calling Pipeline
//
// Pipeline stage definitions for exome variant calling 
// pipeline. See pipeline.groovy for more information.
//
// Author: Simon Sadedin, MCRI
//         Members of Melbourne Genomics
//
// Copyright Melbourne Genomics Health Alliance members. All rights reserved.
//
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
/// 
////////////////////////////////////////////////////////

ENABLE_CADD=true

set_target_info = {

    doc "Validate and set information about the target region to be processed"

    branch.batch = batch
    branch.target_name = branch.name
    branch.target_bed_file = "../design/${target_name}.bed"
    branch.target_gene_file = "../design/${target_name}.genes.txt"
    branch.target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample

    println "Checking for target bed file : $target_bed_file"

    if(!file(target_bed_file).exists() || !file(target_gene_file).exists()) {
        exec """
                mkdir -p ../design;
                cp $BASE/designs/flagships/${target_name}.bed $target_bed_file; 
                cp $BASE/designs/flagships/${target_name}.genes.txt $target_gene_file

                if [ -e $BASE/designs/flagships/${target_name}.pgx.vcf ];
                then
                    cp $BASE/designs/flagships/${target_name}.pgx.vcf ../design;
                fi
        """
    }

    if(!file(target_gene_file).exists())
        fail("Target gene category file $gene_file could not be located for processing sample $branch.name")

    if(!file(target_bed_file).exists())
        fail("Target bed file $target_bed_file could not be located for processing sample $branch.name")

    println "Target $target_name is processing samples $target_samples"
}

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    if(sample_info[sample].target != target_name) {
        // This is expected because every file is processed for every target/flagship
        succeed "No files to process for sample $sample, target $target_name"
    }

    // Patient specific variants are not supported yet
    // If they are provided, we should not process the patient at all
    check {
        if(sample_info[sample].variantsFile?.trim()) 
                exec "false" // force failure
    } otherwise { 
        succeed """
             Study $sample is configured with a sample specific variant file. The pipeline currently does 
             not support sample specific variants. Please remove the variant file from the configuration
             to allow processing.
         """.trim().stripIndent() to channel: gmail, subject: "Invalid configuration for Study $sample"
    }

    def files = sample_info[sample].files.fastq

    println "Processing input files ${files} for target region ${target_bed_file}"
    forward files
}

check_tools = {
    doc """
        Checks for presence of optional tools and sets appropriate pipeline variables
        to enable or disable corresponding pipeline features
        """

    def CADD_DATA = "$BASE/tools/annovar/humandb/hg19_cadd.txt"
    if(!file(CADD_DATA).exists()) {
        msg "Unable to locate data for CADD annotations: CADD scores will not be included"
        ENABLE_CADD = false
    }

    produce("revision.txt") {
        exec """
            git describe --always > $output.txt || true
        """
    }
}

check_sample_info = {

    doc "Validate basic sample information is correct"

    def missingSummary = []
    for(sample in samples) {

        // Check that FASTQ files begin with the sample name followed by underscore
        def files = sample_info[sample].files.fastq
        if(files.any { !file(it).name.startsWith(sample+"_")}) {
            files.each { println "FASTQ: $it | sample=$sample" }
            fail report('templates/invalid_input.html') to channel: gmail, 
                                                           subject: "FASTQ files for sample $sample have invalid file name format", 
                                                           message: "Files $files do not start with the sample name $sample" 
        }

        // Check that all the files specified for the sample exist
        def missingFiles = files.grep { !file(it).exists() }
        if(missingFiles) 
            missingSummary << """
                    The following files specified for sample $sample could not be found:\n\n${missingFiles*.center(120).join('\n')}

                    Please check that the files in your sample file really exist in the data directory.
                 """.stripIndent()

        // Check that file names contain the lane information
        def missingLanes = files.grep { !(it ==~ ".*_L[0-9]*_.*") }
        if(missingLanes) 
            missingSummary << """
                    The following files specified for sample $sample do not contain lane information:\n\n${missingLanes*.center(120).join('\n')}

                    FASTQ file names are required to contain lane identifiers such as L001, L1 or similar. 
                    Please check your input FASTQ and rename it if necessary.
            """

        // Check that file names contain the read number information
        def missingRP = files.grep { !(it ==~ ".*_R[0-9][_.].*fastq.gz\$") }
        if(missingRP) 
            missingSummary << """
                    The following files for sample $sample do not contain the read number in the expected format:\n\n${missingRP*.center(120).join('\n')}

                    FASTQ file names are required to contain the number of the read from the read pair (1 or 2) 
                    in the form '_R1_' or '_R1.'. Please check your input FASTQ and rename it if necessary.
            """
    }

    if(missingSummary) {
        fail missingSummary.join("\n" + ("-" * 120) + "\n")
    }
}

create_combined_target = {

    // Construct the region for variant calling from 
    //
    //   a) all of the disease cohort BED files
    //   b) the EXOME target regions
    //
    // This way we avoid calling variants over the entire genome, but still
    // include everything of interest
    String diseaseBeds = DISEASE_COHORTS.split(",")*.trim().collect { it+".bed"}.join(",")

    output.dir = "../design"

    produce("combined_target_regions.bed") {
        exec """
            cat $BASE/designs/flagships/{$diseaseBeds} $EXOME_TARGET | 
                cut -f 1,2,3 | 
                $BEDTOOLS/bin/bedtools sort | 
                $BEDTOOLS/bin/bedtools merge > $output.bed
        """

        branch.COMBINED_TARGET = output.bed
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "$FASTQC/fastqc -o ${output.dir} $inputs.gz"
    }
}

check_fastqc = {

    doc "Search for any failures in FastQC output and abort further processing if they are found"

    check {

       // NOTE: we remove per-base-sequence content and
       // per-base-gc-content from examination because Nextera
       // appears to contain natural biases that flag QC failures 
       // here.
       exec """
           cat "fastqc/${sample}_*fastqc/summary.txt" |
               grep -v "Per base sequence content" |
               grep -v "Per base GC content" |
               grep -q 'FAIL' && exit 1

           exit 0
       ""","local"
    } otherwise {
        succeed report('templates/fastqc_failure.html') to channel: gmail, 
                                                        subject: "Sample $sample has failed FastQC Check", 
                                                        file: input.zip
    }

    check("FASTQ Format") {
        exec """
            awk -F'\\t' '/Illumina/ { match(\$2, /[0-9.]+/ , result); exit(result[0]<1.7) }' fastqc/${sample}_*_fastqc/fastqc_data.txt
        ""","local"
    } otherwise {
        println "=" * 100
        println "Sample $sample is encoded using a quality encoding incompatible with this pipeline."
        println "Please convert the data first using maq ill2sanger."
        println "=" * 100

        succeed report('templates/fastqc_failure.html') to channel: gmail, 
                                                        subject: "Sample $sample is encoded with incompatible quality scores (Illumina < 1.7)", 
                                                        file: input.zip
    }
}

align_bwa = {

    doc "Align with bwa mem algorithm."

    output.dir = "align"

    var seed_length : 19

    def lanes = inputs.gz.collect { (it.toString() =~ /_(L[0-9]{1,3})_/)[0][1] }.unique()
    if(lanes.size()!=1) 
        succeed report('templates/invalid_input.html') to channel: gmail, 
                                                          subject: "Invalid input files for sample $sample: Bad lane information",
                                                          message: """Failed to identify a unique lane number from FASTQ files: ${inputs.gz}. 
                                                                        Please check the format of the input file names""".stripIndent()
        
    branch.lane = lanes[0]

    def outputFile = sample + "_" + Hash.sha1(inputs.gz*.toString().join(",")) + "_" + lane + ".bam"
    produce(outputFile) {
        //    Note: the results are filtered with flag 0x100 because bwa mem includes multiple 
        //    secondary alignments for each read, which upsets downstream tools such as 
        //    GATK and Picard.
        exec """
                set -o pipefail

                $BWA mem -M -t $threads -k $seed_length 
                         -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample_info[sample].library}\\tSM:${sample}"  
                         $REF $input1.gz $input2.gz | 
                         $SAMTOOLS/samtools view -F 0x100 -bSu - | $SAMTOOLS/samtools sort - $output.prefix
        ""","bwamem"
    }
}

merge_pgx = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"

    if(!file("../design/${target_name}.pgx.vcf").exists()) {
        forward input.recal.vcf.toString() // workaround for Bpipe bug
        return
    }

    msg "Merging vcf files: " + inputs.vcf
    exec """
            java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant $input.recal.vcf
            --variant $input.pgx.vcf
            --out $output.vcf
         """
}

merge_vcf = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"

    produce(target_name + ".merge.vcf") {
        msg "Merging vcf files: " + inputs.vcf
        exec """
                java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
                -T CombineVariants
                -R $REF
                ${inputs.vcf.withFlag("--variant")}
                --out $output.vcf
             """
    }
}

merge_bams = {
    doc """
        Merge the BAM files from multiple lanes together.
        <p>
        When there is only one BAM file the merge is skipped, and a 
        symbolic link created instead.
        """

    output.dir="align"
    produce(sample + ".merge.bam") {

        // If there is only 1 bam file, then there is no need to merge
        if(inputs.bam.size()==1)  {
                msg "Skipping merge of $inputs.bam because there is only one file"
                // This use of symbolic links may be questionable
                // However if the ordinary case involves only one
                // bam file then there may be some significant savings
                // from doing this.
                exec "ln -s ${file(input.bam).name} $output.bam; ln -s ${file(input.bam).name}.bai ${output.bam}.bai;"
        }
        else {
            msg "Merging $inputs.bam size=${inputs.bam.size()}"
            exec """
                java -Xmx2g -jar $PICARD_HOME/lib/MergeSamFiles.jar
                    ${inputs.bam.withFlag("INPUT=")}
                    VALIDATION_STRINGENCY=LENIENT
                    ASSUME_SORTED=true
                    CREATE_INDEX=true
                    OUTPUT=$output.bam
             """, "merge"
        }
    }
}

index_bam = {

    doc "Create an index for a BAM file"

    // A bit of a hack to ensure the index appears in the
    // same directory as the input bam, no matter where it is
    // nb: fixed in new version of Bpipe
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "$SAMTOOLS/samtools index $input.bam"
    }
    forward input
}

realignIntervals = {
    doc "Discover candidate regions for realignment in an alignment with GATK"
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T RealignerTargetCreator 
            -R $REF 
            -I $input.bam 
            --known $GOLD_STANDARD_INDELS 
            -o $output.intervals
    """, "realign_target_creator"
}

realign = {
    doc "Apply GATK local realignment to specified intervals in an alignment"
    output.dir="align"
    exec """
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T IndelRealigner 
             -R $REF 
             -I $input.bam 
             -targetIntervals $input.intervals 
             -o $output.bam
    ""","local_realign"
}

dedup = {
    doc "Remove PCR duplicates from reads"
    output.dir="align"
    exec """
        java -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $PICARD_HOME/lib/MarkDuplicates.jar
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$output.metrics
             OUTPUT=$output.bam
    """
}

recal_count = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    INDEL_QUALS=""
    // To use lite version of GATK uncomment below
    // INDEL_QUALS="--disable_indel_quals"
    exec """
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T BaseRecalibrator 
             -I $input.bam 
             -R $REF 
             --knownSites $DBSNP $INDEL_QUALS
             -l INFO 
             -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate 
             -o $output.counts
    """, "recalibrate_bam"
}

recal = {
    doc "Apply recalibration quality adjustments so that quality scores match actual observed error rates"
    output.dir="align"
    exec """
          java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
               -T PrintReads 
               -I $input.bam 
               -BQSR $input.counts 
               -R $REF 
               -l INFO 
               -o $output.bam
        """, "recalibrate_bam"
}

call_variants = {
    doc "Call SNPs/SNVs using GATK Unified Genotyper"
    output.dir="variants"

    // Default values of confidence thresholds
    // come from the Broad web site. However
    // these may be higher than suitable in our context
    var call_conf:5.0, 
        emit_conf:5.0

    transform("bam","bam") to("metrics","vcf") {
        exec """
                java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
                   -R $REF 
                   -I $input.bam 
                   -nt 4
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   -dcov 1600 
                   -l INFO 
                   -L $COMBINED_TARGET
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -glm BOTH
                   -metrics $output.metrics
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

call_pgx = {
    doc "Call Pharmacogenomic variants using GATK Unified Genotyper"
    output.dir="variants"

    var call_conf:5.0, 
        emit_conf:5.0

    if(!file("../design/${target_name}.pgx.vcf").exists())
        return

    transform("bam","bam") to("metrics","pgx.vcf") {
        exec """
                java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
                   -R $REF 
                   -I $input.bam 
                   -nt 4
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   --output_mode EMIT_ALL_SITES
                   -dcov 1600 
                   -l INFO 
                   -L ../design/${target_name}.pgx.vcf
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -glm BOTH
                   -metrics $output.metrics
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

@filter("filter")
filter_variants = {
    doc "Select only variants in the genomic regions defined for the $target_name target"
    output.dir="variants"

    def pgx_flag = ""
    if(file("../design/${target_name}.pgx.vcf").exists()) {
        pgx_flag = "-L ../design/${target_name}.pgx.vcf"
    }

    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
             -R $REF
             -T SelectVariants 
             --variant $input.vcf 
             -L $target_bed_file $pgx_flag
             -o $output.vcf 
    """
}

@filter("vep")
annotate_vep = {
    doc "Annotate variants using VEP to add Ensemble annotations"
    output.dir="variants"
    exec """
        PERL5LIB="$CONDEL"
        perl $VEP/variant_effect_predictor.pl --cache --dir $TOOLS/vep/vep_cache 
            -i $input.vcf 
            --vcf -o $output.vcf 
            -species human 
            --canonical --per_gene --protein 
            --sift=b --polyphen=b
            --symbol hgnc --force_overwrite --hgvs  --maf_1kg --maf_esp --pubmed
            --plugin Condel,$CONDEL/config,s
    """, "vep"
}

calc_coverage_stats = {
    doc "Calculate coverage across a target region using Bedtools"
    output.dir="qc"
    transform("bam","bam") to(file(target_bed_file).name+".cov.txt","ontarget.txt") {
        exec """
          $BEDTOOLS/bin/coverageBed -d  -abam $input.bam -b $target_bed_file > $output.txt

          $SAMTOOLS/samtools view -L $COMBINED_TARGET $input.bam | wc | awk '{ print \$1 }' > $output2.txt
        """
    }
}

check_coverage = {

    output.dir = "qc"

    def medianCov
    transform("cov.txt") to("cov.stats.median", "cov.stats.csv") {

        R {"""
            bam.cov = read.table("$input.cov.txt", col.names=c("chr","start","end", "gene", "offset", "cov"))
            meds = aggregate(bam.cov$cov, list(bam.cov$gene), median)
            write.csv(data.frame(Gene=meds[,1],MedianCov=meds$x), "$output.csv", quote=F, row.names=F)
            writeLines(as.character(median(bam.cov$cov)), "$output.median")
        """}

        medianCov = Math.round(file(output.median).text.toFloat())
    }

    check {
        exec "[ $medianCov -ge $MEDIAN_COVERAGE_THRESHOLD ]"
    } otherwise {
        // It may seem odd to call this a success, but what we mean by it is that
        // Bpipe should not fail the whole pipeline, merely this branch of it
        succeed report('templates/sample_failure.html') to channel: gmail, 
                                                           median: medianCov, 
                                                           file:output.csv, 
                                                           subject:"Sample $sample has failed with insufficient median coverage ($medianCov)"
    }
}

check_karyotype = {

    doc "Compare the inferred sex of the sample to the inferred karyotype from the sequencing data"

    def karyotype_file = "results/" + sample + '.summary.karyotype.tsv'
    check {
        exec """
            [ `grep '^Sex' $karyotype_file | cut -f 2` == `grep 'Inferred Sex' $karyotype_file | cut -f 2` ]
        """
    } otherwise {
        // It may seem odd to call this a success, but what we mean by it is that
        // Bpipe should not fail the whole pipeline, merely this branch of it
        succeed report('templates/sample_failure.html') to channel: gmail, 
                                                           median: medianCov, 
                                                           file: karyotype_file,
                                                           subject:"Sample $sample has a different sex than inferred from sequencing data"
     }
}

augment_condel = {

    doc "Extract Condel scores from VEP annotated VCF files and add them to Annovar formatted CSV output"

    output.dir="variants"
    from("exome_summary.csv") filter("con") {
        exec """
            JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY 
                -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar 
                $SCRIPTS/merge_condel.groovy
                    -i $input.vcf
                    -a $input.csv
                    -o $output.csv
        """
    }
}

augment_cadd = {

    doc "Add CADD annotations to an existing file of Annovar annotations"

    if(!ENABLE_CADD) 
        return

    output.dir = "variants"
    from("*.exome_summary.*.csv","*.hg19_cadd_dropped" ) filter("cadd") {
        exec """
            JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY 
                -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar 
                $SCRIPTS/add_cadd_scores.groovy
                    -a $input.csv
                    -c $input.hg19_cadd_dropped
                    -o $output.csv
        """
    }
}

calculate_cadd_scores = {

    doc "Compute CADD scores by running Annovar annotation"

    if(!ENABLE_CADD) 
        return

    output.dir="variants"
    exec """
        $ANNOVAR/annotate_variation.pl
        $input.av 
        $ANNOVAR/../humandb/
        -filter 
        -dbtype cadd 
        -buildver hg19 
        -out $output.hg19_cadd_dropped.prefix
    """
}

index_vcf = {
    output.dir="variants"
    transform("vcf","vcf") to("sort.vcf","sort.vcf.idx") {
        exec """
            $IGVTOOLS/igvtools sort $input.vcf $output.vcf 

            $IGVTOOLS/igvtools index $output.vcf
        """
    }
}

vcf_to_excel = {

    doc "Convert a VCF output file to Excel format, merging information from Annovar"

    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    check {
        exec "ls results/${target_name}.qc.xlsx > /dev/null 2>&1"
    } otherwise { 
        succeed "No samples succeeded for target $target_name" 
    }

    def pgx_flag = ""
    if(file("../design/${target_name}.pgx.vcf").exists()) {
        pgx_flag = "-pgx ../design/${target_name}.pgx.vcf"
    }

    output.dir="results"

    def all_outputs = [target_name + ".xlsx"] + target_samples.collect { it + ".annovarx.csv" }
    from("*.exome_summary.*.csv", "*.vcf") produce(all_outputs) {
        exec """
            echo "Creating $outputs.csv"

            JAVA_OPTS="-Xmx2g -Djava.awt.headless=true" $GROOVY 
                -cp $SCRIPTS:$GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar:$TOOLS/sqlite/sqlitejdbc-v056.jar $SCRIPTS/vcf_to_excel.annovar.groovy 
                -s '${target_samples.join(",")}'
                ${inputs.csv.withFlag("-a")}
                ${inputs.vcf.withFlag("-vcf")}
                -x "synonymous SNV"
                -db $VARIANT_DB
                -o $output.xlsx
                -si $sample_metadata_file
                -gc $target_gene_file ${pgx_flag}
                -annox $output.dir
                ${inputs.bam.withFlag("-bam")}
        """
    }
}

plot_coverage = {
    doc "Create plots showing coverage distributions for alignment"
    from("cov.txt") {
        transform("cum.png") {
            msg "Plotting coverage for sample $sample_name"
            R {"""
                sample = '$input.txt'
                name = '$sample'
                samplecov = read.table(file=pipe(paste("awk '{ print $NF }' ",sample, sep='')))
                tf = table(factor(samplecov$V1, levels=0:max(samplecov)))
                cs = cumsum(tf)
                png("$output.png")
                plot(1-cs/max(cs), ylim=c(0,1), type="l", xaxs="i", xlim=c(0,800),
                     xaxt="n",
                     main=paste("Cum. Coverage Distribution for ", name), xlab='Coverage Depth', ylab='Frequency')
                axis(side=1, at=seq(0,800,25))
                dev.off()
            """} 
        }
    }
}

gatk_depth_of_coverage = {

    doc "Calculate statistics about depth of coverage for an alignment using GATK"

    output.dir = "qc"
    transform("bam") to("."+target_name + ".cov.sample_cumulative_coverage_proportions", 
                         "."+target_name + ".cov.sample_interval_statistics") { 
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
               -R $REF
               -T DepthOfCoverage 
               -o $output.sample_cumulative_coverage_proportions.prefix
               --omitDepthOutputAtEachBase
               -I $input.bam
               -ct 1 -ct 10 -ct 20 -ct 50 -ct 100
               -L $target_bed_file
        """
    }
}

qc_excel_report = {

    doc "Create an excel file containing a summary of QC data for all the samples for a given target region"

    var LOW_COVERAGE_THRESHOLD : 15,
        LOW_COVERAGE_WIDTH : 1

    output.dir="results"

    def samples = sample_info.grep { it.value.target == target_name }.collect { it.value.sample }
    from("*.cov.txt", "*.dedup.metrics") produce(target_name + ".qc.xlsx") {
            exec """
                JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/qc_excel_report.groovy 
                    -s ${target_samples.join(",")} 
                    -t $LOW_COVERAGE_THRESHOLD
                    -w $LOW_COVERAGE_WIDTH
                    -o $output.xlsx
                    $inputs.sample_cumulative_coverage_proportions  
                    $inputs.sample_interval_statistics 
                    $inputs.metrics 
                    $inputs.txt
            ""","qc_excel_report"
    }
}

@filter("sig")
annotate_significance = {
    doc "Add clinical significance category annotations as defined by Melbourne Genomics"
    output.dir="variants"
    from("con.csv") {
        exec """
            python $SCRIPTS/annotate_significance.py -a $input.csv > $output.csv
        """
    }
}

annovar_summarize_refgene = {
    doc "Annotate variants using Annovar"
    output.dir="variants"
    transform("vcf") to(["av","av.refgene.exome_summary.csv","av.refgene.exonic_variant_function","av.refgene.genome_summary.csv"]) {
        exec """
            $ANNOVAR/convert2annovar.pl $input.vcf -format vcf4 > $output.av

            $ANNOVAR/summarize_annovar.pl 
                --genetype refgene 
                --verdbsnp 138  
                --outfile ${output.av}.refgene 
                --buildver hg19  $output.av $ANNOVAR/../humandb/
        """
    }
}

add_to_database = {
    doc "Add discovered variants to a database to enable annotation of number of observations of the variant"
    output.dir="variants"
    uses(variantdb:1) {
        exec """

            echo "====> Adding variants for flaship $target_name to database"

            JAVA_OPTS="-Xmx24g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar:$TOOLS/sqlite/sqlitejdbc-v056.jar $SCRIPTS/vcf_to_db.groovy 
                   -v $input.recal.vcf 
                   -a $input.csv 
                   -db $VARIANT_DB 
                   -cohort $target_name
                   -b "$batch"

            echo "<==== Finished adding variants for flaship $target_name to database"

            echo "Variants from $input.vcf were added to database $VARIANT_DB on ${new Date()}" > $output.txt
        """, "add_to_database"
    }
}

reorder = {
    filter('reorder') {
        exec """
            java -Xmx2g -jar $PICARD_HOME/lib/ReorderSam.jar
                TMP_DIR=$TMP_DIR
                I=$input.bam
                O=$output.bam
                VALIDATION_STRINGENCY=LENIENT
                REFERENCE=$REF
            """ 
    }
}

summary_pdf = {

    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"

    produce("${sample}.summary.pdf","${sample}.summary.karyotype.tsv") {
        exec """
             JAVA_OPTS="-Xmx3g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar $SCRIPTS/qc_pdf.groovy 
                -cov $input.cov.txt
                -ontarget $input.ontarget.txt
                -metrics $input.metrics
                -study $sample 
                -meta $sample_metadata_file
                -threshold 20 
                -classes GOOD:95:GREEN,PASS:80:ORANGE,FAIL:0:RED 
                -exome $EXOME_TARGET
                -gc $target_gene_file ${pgx_flag}
                -bam $input.bam
                -o $output.pdf
        """

        branch.karyotype = output.tsv

        send text {"Sequencing Results for Study $sample"} to channel: gmail, file: output.pdf
    }
}

sample_similarity_report = {

    doc "Create a report indicating the difference in count of variants for each combination of samples"

    output.dir = "qc"

    produce("similarity_report.txt") {
        exec """
            (
                for i in ${samples.join(' ')} ;
                do  
                    for j in ${samples.join(' ')} ;
                    do 
                        if [ "$i" != "$j" ]; then
                            printf "$i\\t$j\\t"; 

                            JAVA_OPTS=-Xmx1g $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar $SCRIPTS/vcf_diff.groovy 
                                    -vcf variants/\$i.merge.dedup.realign.recal.vcf 
                                    -vcf variants/\$j.merge.dedup.realign.recal.vcf | wc | awk '{ print \$1 }'; echo ; 
                        fi
                    done; 
                done
            ) > $output.txt
        """
    }
}

provenance_report = {
    branch.sample = branch.name
    output.dir = "results"
    produce(sample + ".provenance.pdf") {
       send report("scripts/provenance_report.groovy") to file: output.pdf
    }
}

annovar_to_lovd = {
    branch.sample = branch.name
    output.dir="results/lovd"
    produce(sample +"_LOVD") {
        exec """
            python scripts/annovar2LOVD.py --csv $input.annovarx.csv --meta $sample_metadata_file --outdir results/lovd
        """
    }
}

