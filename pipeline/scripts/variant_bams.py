#!/usr/bin/env python2.7
####################################################################################
#
# Melbourne Genomics Variant BAM generator
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
# DISTRIBUTION:
#
# This source code should not be distributed to a third party without prior
# approval of the Melbourne Genomics Health Alliance steering committee (via
# Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
#
####################################################################################
#
# Purpose:
# * Generate bam files from a list of variants
####################################################################################

import os
import csv
from subprocess import call
from argparse import (ArgumentParser, FileType, ArgumentDefaultsHelpFormatter)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Produce a bam file of reads overlapping a variant. By default, reads overlapping the region 100bp upstream and downstream of the variant are included.',
    formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--bam', type=str, required=True,
        help='Input bam file.')
    parser.add_argument(
        '--csv', type=str, required=True,
        help='Variants in csv format. Must have at least the columns AAChange, Chr, Start, End. Sample assumed to be the start of this file name (before the ".").')
    parser.add_argument(
        '--outdir', type=str, default='.',
        help='Output directory for bams. Directory must already exist.')
    parser.add_argument(
        '--log', type=str,
        help='Optionally output the number of bams created to a log file with this name.')
    parser.add_argument(
        '--upstream', type=int, default=100,
        help="Include reads this many bp upstream (5') of the variant.")
    parser.add_argument(
        '--downstream', type=int, default=100,
        help="Include reads this many bp downstream (3') of the variant.")
    parser.add_argument(
        '--samtoolsdir', type=str,
        help="Directory in which samtools is installed. Default: samtools assumed to be in PATH.")

    return parser.parse_args() 


def main():
    # Parse command line arguments
    args = parse_args()
    inbam = args.bam
    variantfile = args.csv
    outdir = args.outdir
    upstream = args.upstream # Default 100 
    downstream = args.downstream # Default 100
    if args.samtoolsdir:
        samtools_exec = args.samtoolsdir + '/samtools'
    else:
        samtools_exec = 'samtools'

    # Assume sample name is the first part of the filename before the "."
    sample = variantfile.split('/')[-1].split('.')[0]

    with open(variantfile) as variantcsv:
        var_count = 0 
        for line in csv.DictReader(variantcsv):
            NM = line['AAChange'].split(':')[0]

            chr = line['Chr']
            start = int(line['Start'])
            end = int(line['End'])

            outbam = '{0}-{1}-{2}-{3}-{4}-IGV.bam'.format(sample, NM, chr, start, end)
            region = '{0}:{1}-{2}'.format(chr, start - upstream, end + downstream)

            # create new bam containing reads in the given region
            call([samtools_exec, 'view', '-b', '-o', outdir+'/'+outbam, inbam, region])
            # index the bam
            call([samtools_exec, 'index', outdir+'/'+outbam])
            
            var_count += 1
            if var_count % 100 == 0:
              print "{0} variants processed...".format( var_count )

    # If required, produce a log file
    # This is useful to trick bpipe into tracking the output of this script, 
    # since an unknown number of bam files are produced (sometimes zero)
    if args.log:
        with open(args.log, 'w') as logfile:
            logfile.write('{0} bams for sample {1} written to {2}'.format(var_count, sample, outdir))

if __name__ == '__main__':
    main()