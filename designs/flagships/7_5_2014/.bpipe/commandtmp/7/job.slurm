#!/bin/bash
#SBATCH --job-name=annotate
#SBATCH --account VR0002
#SBATCH --mem=8192
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH -p main

set -o errexit

 (echo "Annotating genes ...";          JAVA_OPTS="-Xmx1g" /vlsci/VR0320/shared/production/tools/groovy/2.1.1/bin/groovy /vlsci/VR0320/shared/production/pipeline/scripts/annotate_genes.groovy -r /vlsci/VR0320/shared/production/tools/annovar/humandb/hg19_refGene.txt MelbourneGenomics-CCDS-CRC-20140507.flatten.bed > MelbourneGenomics-CCDS-CRC-20140507.flatten.annotate.bed) > .bpipe/commandtmp/7/7.out 2>  .bpipe/commandtmp/7/7.err
