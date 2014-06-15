#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account VR0002
#SBATCH --mem=8192
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH -p main

set -o errexit

 (cp MelbourneGenomics-CCDS-EPIL-20140507.flatten.annotate.sort.bed EPIL.bed) > .bpipe/commandtmp/17/17.out 2>  .bpipe/commandtmp/17/17.err
