#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account VR0002
#SBATCH --mem=8192
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH -p main

set -o errexit

 (cp MelbourneGenomics-CCDS-CMT-20140507.flatten.annotate.sort.bed CMT.bed) > .bpipe/commandtmp/18/18.out 2>  .bpipe/commandtmp/18/18.err
