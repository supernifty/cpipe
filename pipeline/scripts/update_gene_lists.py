####################################################################################
#
# Melbourne Genomics Pipeline Annotation Script
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
#   Given a gene list and reference bed file, generate a bed file with just those genes
# Usage:
#   update_gene_lists --source dir --target dir
####################################################################################

import argparse
import datetime
import glob
import os
import re
import sys

CATEGORY = '1'

def write_log( log, msg ):
  log.write( '%s: %s\n' % ( datetime.datetime.now().strftime( '%y%m%d-%H%M%S' ), msg ) )

def update_gene_lists( source_dir, target_dir, log ):
  for filename in glob.glob( os.path.join( source_dir, '*.add.genes.txt' ) ):
    cohort = os.path.basename( filename ).split( '.' )[0]
    # find corresponding flagship
    target = os.path.join( target_dir, '%s.genes.txt' % cohort )
    if os.path.isfile( target ):
      # read existing genes and categories
      genes = {}
      for line in open( target, 'r' ):
        if line.startswith('#'):
          continue
        fields = line.strip().split('\t')
        genes[fields[0].upper()] = fields[1]
      # read new genes
      added = set()
      candidates = set()
      for line in open( filename, 'r' ):
        if line.startswith('#'):
          continue
        gene = line.strip().upper()
        candidates.add( gene )
        if gene not in genes:
          added.add(gene)
          genes[gene] = CATEGORY
        
      # write out additional
      if len(added) > 0:
        with open( target, 'w' ) as fh:
          fh.write( '#version %s\n' % datetime.datetime.now().strftime( '%y%m%d' ) )
          fh.write( '#notes %i gene(s) added: %s\n' % ( len(added), ','.join( sorted( list(added) ) ) ) )
          for gene in sorted(genes.keys()):
            fh.write( '%s\t%s\n' % ( gene, genes[gene] ) )
        write_log( log, '%s: %i gene(s) added from %i candidate(s): %s' % ( cohort, len(added), len(candidates), ','.join( sorted( list(added) ) ) ) )
      else:
        write_log( log, '%s: no changes from %i candidate(s)' % ( cohort, len(candidates) ) )
    else:
      write_log( log, 'ERROR: target gene list %s does not exist' % target )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate bed files')
  parser.add_argument('--source', required=True, help='source of extra genes') # input
  parser.add_argument('--target', required=True, help='target containing gene files to update') # input
  parser.add_argument('--log', required=False, help='write changes to this file') # input
  args = parser.parse_args()
  if args.log:
    log = open( args.log, 'a+' )
  else:
    log = sys.stderr

  update_gene_lists( args.source, args.target, log )
