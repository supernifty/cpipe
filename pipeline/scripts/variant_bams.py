#!/usr/bin/env python2.7
import os
import csv
from subprocess import call
from argparse import (ArgumentParser, FileType)

# module load samtools-intel/0.1.19

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Produce a bam file of reads overlaping a variant')
    parser.add_argument(
        '--bam', type=str, required=True,
        help='Input bam file')
    parser.add_argument(
        '--csv', type=str, required=True,
        help='Variants in csv format. Must have at least the columns AAChange, Chr, Start, End. Sample assumed to be start of this file name (before the .)')

    return parser.parse_args() 

def main():
    # Parse command line arguments
    args = parse_args()
    inbam = args.bam
    variantfile = args.csv

# # Specific sample/variant test:
# # 020432001
# # Gene: 'SLC19A3'
# baseDir = '/vlsci/VR0320/shared/hdashnow/MelbGenomicsPipeline/batches/007_development'
# variantfile = baseDir + '/analysis/results/510000101.annovarx.csv'
# inbam = baseDir + '/analysis/align/510000101.merge.dedup.realign.recal.bam'
# 
# # 020432701
# # Gene: ATP1A3
# baseDir = '/vlsci/VR0320/shared/production/batches/008'
# variantfile = baseDir + '/results/020432701.annovarx.csv'
# inbam = baseDir + '/analysis/align/020432701.merge.dedup.realign.recal.bam'

    # These should be set at the command line
    upstream = 100
    downstream = 100
    sample = variantfile.split('/')[-1].split('.')[0]

    with open(variantfile) as variantcsv:
        for line in csv.DictReader(variantcsv):
            NM = line['AAChange'].split(':')[0]

#            # For testing
#            if line['Gene'] != 'ATP1A3':
#                continue

            chr = line['Chr']
            start = int(line['Start'])
            end = int(line['End'])
    #        print NM, chr, start, end 

            outbam = '{0}-{1}-{2}-{3}-{4}-IGV.bam'.format(sample, NM, chr, start, end)
            region = '{0}:{1}-{2}'.format(chr, start - upstream, end + downstream)

            call(['samtools', 'view', '-b', '-o', outbam, inbam, region])


if __name__ == '__main__':
    main()
