// vim: set ts=4:sw=4:expandtab:cindent
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
// License: TODO
// 
////////////////////////////////////////////////////////

set_target_info = {

    doc "Validate and set information about the target region to be processed"

    branch.target_name = branch.name
    branch.target_bed_file = "../design/${target_name}.bed"
    branch.target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample

    if(!new File(target_bed_file).exists())
        fail("Target bed file $target_bed_file could not be located for processing sample $branch.name")

    println "Target $target_name is processing samples $target_samples"
}

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name

    if(sample_info[sample].target != target_name) {
        succeed "No files to process for sample $sample, target $target_name"
    }
    else {
        def files = sample_info[sample].files
        println "Processing input files ${files} for target region ${target_bed_file}"
        forward files
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "fastqc -o ${output.dir} $inputs.gz"
    }
}

align_bwa = {
    doc """
        Align with bwa mem algorithm.
    """

    //    Note: the results are filtered with flag 0x100 because bwa mem includes multiple 
    //    secondary alignments for each read, which upsets downstream tools such as 
    //    GATK and Picard.
    output.dir = "align"

    var seed_length : 19

    // TODO: replace with real platform unit, ID and lane value based on real 
    // meta data file. For now these are faked
    produce(sample + ".bam") {
        exec """
                set -o pipefail

                $BWA mem -M -t $threads -k $seed_length 
                         -R "@RG\\tID:1\\tPL:illumina\\tPU:1\\tLB:1\\tSM:${sample}"  
                         $REF $input1.gz $input2.gz | 
                         samtools view -F 0x100 -bSu - | samtools sort - $output.prefix
        ""","bwamem"
    }
}

merge_vcf = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"
    produce(target_name + ".merge.vcf") {
        msg "Merging vcf files: " + inputs.vcf
        exec """
                java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
                -T CombineVariants
                -R $HGFA
                -L $target_bed_file
                ${inputs.vcf.withFlag("--variant")}
                --out $output.vcf
             """
    }
}

index_bam = {
    // A bit of a hack to ensure the index appears in the
    // same directory as the input bam, no matter where it is
    // nb: fixed in new version of Bpipe
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

realignIntervals = {
    // Hard-coded to take 2 known indels files right now
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I $input.bam --known $GOLD_STANDARD_INDELS -o $output.intervals
    """
}

realign = {
    output.dir="align"
    exec """
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $input.bam -targetIntervals $input.intervals -o $output.bam
    ""","local_realign"
}

dedup = {
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
    // TODO: use non-lite version of GATK, remove disable_indel_quals flag
    exec """
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T BaseRecalibrator 
             -I $input.bam 
             -R $REF 
             --knownSites $DBSNP 
             --disable_indel_quals
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
    var call_conf:50.0, 
        emit_conf:10.0

    transform("bam","bam") to("metrics","vcf") {
        exec """
                java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
                   -R $REF 
                   -I $input.bam 
                   -nt 4
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   -dcov 1600 
                   -L $EXOME_TARGET
                   -l INFO 
                   -A AlleleBalance -A DepthOfCoverage -A FisherStrand 
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
    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
             -R $REF
             -T SelectVariants 
             --variant $input.vcf 
             -L $target_bed_file
             -o $output.vcf 
    """
}

@filter("vep")
annotate = {
    doc "Annotate variants using VEP to add Ensemble annotations"
    output.dir="variants"
    exec """
        perl $VEP/variant_effect_predictor.pl --cache --dir $VEP/vep_cache 
            -i $input.vcf 
            --vcf -o $output.vcf 
            -species human 
            --canonical --per_gene --protein 
            --sift=b --polyphen=b
            --symbol hgnc --force_overwrite --hgvs  --maf_1kg --maf_esp --pubmed
    """
}

calc_coverage_stats = {
    doc "Calculate coverage across a target region"
    output.dir="qc"
    from(target_bed_file) {
        transform("bam") to(file(input.bed).name+".cov.txt") {
            exec """
              coverageBed -d  -abam $input.bam -b $input.bed > $output.txt
            """
        }
    }
}

sort_vcf = {
    output.dir="variants"
    filter("sort") {
        // NOTE: bedtools sorts in order chr1,chr10,... which is less intuitive
        // and incompatible with reference
        // exec "$BEDTOOLS/bedtools sort -header -i $input.vcf > $output.vcf"
        exec "$IGVTOOLS/igvtools sort $input.vcf $output.vcf"
    }
}

@transform("xlsx")
vcf_to_excel_vep = {
    doc "Convert a VCF output file to Excel format, including annotations from VEP"
    exec """
        JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" groovy 
            -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar 
            $BASE/pipeline/scripts/vcf_to_excel_vep.groovy 
                -i $input.vcf
                -o $output.xlsx
                -t ${new File("..").absoluteFile.name}
    """
}

@transform("xlsx")
vcf_to_excel = {
    doc "Convert a VCF output file to Excel format, merging information from Annovar"

    from("*.exome_summary.csv", "*.vcf") produce(target_name + ".xlsx") {
        exec """
            JAVA_OPTS="-Xmx2g -Djava.awt.headless=true" groovy 
                -cp $EXCEL/excel.jar $BASE/pipeline/scripts/vcf_to_excel.annovar.groovy 
                -s '${target_samples.join(",")}'
                -a $input.csv 
                -i $input.vcf
                -x "synonymous SNV"
                -o $output.xlsx
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

    output.dir = "qc"

    //var target_name : "all"

    transform(".bam") to("."+target_name + ".cov.sample_cumulative_coverage_counts") {
        exec """
            java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
               -R $REF
               -T DepthOfCoverage 
               -o $output.sample_cumulative_coverage_counts.prefix
               -I $input.bam
               -ct 1 -ct 10 -ct 20 -ct 50 -ct 100
               -L $target_bed_file
        """
    }
}

qc_excel_report = {

    doc "Create an excel file containing a summary of QC data for all the samples for a given target region"

    def samples = sample_info.grep { it.value.target == target_name }.collect { it.value.sample }

    from(target_name+".bed.cov.txt") {
        exec """
            JAVA_OPTS=-Xmx4g groovy -cp $EXCEL/excel.jar $BASE/pipeline/scripts/qc_excel_report.groovy 
                -s ${samples.keySet().join(",")} 
                $inputs.sample_cumulative_coverage_proportions  
                $inputs.sample_interval_statistics 
                $inputs.metrics 
                $inputs.txt
        """
    }
}


annovar_summarize_refgene = {
    doc "Annotate variants using Annovar"
    output.dir="variants"
    transform("vcf") to(["av","av.refgene.exome_summary.csv","av.refgene.exonic_variant_function","av.refgene.genome_summary.csv"]) {
        exec """
                $ANNOVAR/convert2annovar.pl $input -format vcf4 > $output.av

                $ANNOVAR/summarize_annovar.pl 
                    --genetype refgene 
                    --verdbsnp 138  
                    --outfile ${output.av}.refgene 
                    --buildver hg19  $output.av $ANNOVAR/../humandb/
        """
    }
}
