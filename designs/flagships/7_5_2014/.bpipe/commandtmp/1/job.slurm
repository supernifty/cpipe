#!/bin/bash
#SBATCH --job-name=flatten
#SBATCH --account VR0002
#SBATCH --mem=8192
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH -p main

set -o errexit

 (echo "Flattening and removing nonstandard chromosomes ...";          sortBed -i raw/MelbourneGenomics-CCDS-CS-20140507.bed | bedtools merge -nms -i - | sed 's/;.*$//' | grep -v '^chr[0-9]_' > MelbourneGenomics-CCDS-CS-20140507.flatten.bed) > .bpipe/commandtmp/1/1.out 2>  .bpipe/commandtmp/1/1.err
