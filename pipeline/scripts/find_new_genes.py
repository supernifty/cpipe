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
#   find_new_genes --reference reference_bed --exclude exclude --target target_dir < sample_file
####################################################################################

import argparse
import collections
import os
import re
import sys

def generate_new_genes( sample_lines, log, reference_genes, excluded_genes ):
  '''
    given samples, make files of the form CS.extra.genes.txt and CS.extra.excluded.genes.txt
  '''
  # get list of reference genes
  reference = set()
  for line in reference_genes:
    fields = line.strip().split('\t')
    if len(fields) > 3:
      reference.add( fields[3].upper() )
  # get list of excluded genes
  excluded = set()
  for line in excluded_genes:
    excluded.add( line.strip().split('\t')[0].upper() )
  log.write( '%i available reference genes, %i excluded genes\n' % ( len(reference), len(excluded) ) )

  # parse sample
  headers = {}
  for idx, title in enumerate(sample_lines[0].split('\t')):
    headers[title] = idx

  if 'Cohort' not in headers:
    log.write( 'ERROR: Cohort not in sample header\n' )     
    return 1
  if 'Prioritised_Genes' not in headers:
    log.write( 'ERROR: Prioritised_Genes not in sample header\n' )     
    return 1
  
  result = {}
  for line in sample_lines[1:]:
    fields = line.split('\t')
    cohort = fields[headers['Cohort']]
    genes = fields[headers['Prioritised_Genes']]
    if cohort != '':
      if cohort not in result:
        result[cohort] = { 'add': set(), 'addonce': set(), 'notfound': set() }
      candidates = re.split( '[:,"]+', genes.upper() )
      for candidate in candidates:
        if candidate in excluded:
          result[cohort]['addonce'].add( candidate )
        elif candidate in reference:
          result[cohort]['add'].add( candidate )
        else:
          if candidate != '' and candidate != '4':
            result[cohort]['notfound'].add( candidate )
  return result

def write_genes( additions, target, log, dummy=False ):
  '''
    writes e.g. CS.add.genes.txt, CS.excluded.genes.txt, CS.notfound.genes.txt
  '''
  for cohort in additions:
    for category in additions[cohort]:
      if dummy:
        for gene in additions[cohort][category]:
          log.write( '%s.%s.%s\n' % (cohort, category, gene ) )
      else:
        with open( os.path.join( target, '%s.%s.genes.txt' % ( cohort, category ) ), 'w' ) as fh:
          for gene in additions[cohort][category]:
            fh.write( '%s\n' % gene )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate bed files')
  parser.add_argument('--reference', required=True, help='reference bed file') # input
  parser.add_argument('--exclude', required=True, help='file containing genes to exclude') # input
  parser.add_argument('--target', required=True, help='target directory')
  args = parser.parse_args()
  samples = sys.stdin.readlines()
  additions = generate_new_genes( samples, sys.stderr, open(args.reference, 'r'), open(args.exclude, 'r') )
  write_genes( additions, args.target, sys.stderr, dummy=False )
