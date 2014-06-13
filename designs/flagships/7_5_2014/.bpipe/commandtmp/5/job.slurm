#!/bin/bash
#SBATCH --job-name=flatten
#SBATCH --account VR0002
#SBATCH --mem=8192
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH -p main

set -o errexit

 (echo "Flattening and removing nonstandard chromosomes ...";          sortBed -i raw/MelbourneGenomics-CCDS-AML-20140507.bed | bedtools merge -nms -i - | sed 's/;.*$//' | grep -v '^chr[0-9]_' > MelbourneGenomics-CCDS-AML-20140507.flatten.bed) > .bpipe/commandtmp/5/5.out 2>  .bpipe/commandtmp/5/5.err
