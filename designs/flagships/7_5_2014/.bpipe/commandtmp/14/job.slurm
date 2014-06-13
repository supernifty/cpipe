#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --account VR0002
#SBATCH --mem=8192
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH -p main

set -o errexit

 (echo "Sorting ...";          /vlsci/VR0320/shared/production/tools/IGVTools/2.3.6/igvtools.lowmem sort MelbourneGenomics-CCDS-CMT-20140507.flatten.annotate.bed MelbourneGenomics-CCDS-CMT-20140507.flatten.annotate.sort.bed) > .bpipe/commandtmp/14/14.out 2>  .bpipe/commandtmp/14/14.err
