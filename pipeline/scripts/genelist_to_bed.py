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
#   Given some gene lists and reference bed file, generate a bed file with just those genes
# Usage:
#   genelist_to_bed genelist... < ref.bed > filtered.bed
#   arguments:
#   --exclude a file containing genes to exclude
####################################################################################

import sys

def filter_bed( genelists, bed_in, bed_out, log, exclude=None ):
  # get the list of proposed genes
  genes = set()
  for arg in genelists:
    for line in open( arg, 'r' ):
      if line.startswith('#'):
        continue
      fields = line.strip().split( '\t' )
      genes.add( fields[0].upper() )
  log.write( '%i candidate genes added\n' % len(genes) )
  
  # get the list of exclusions
  disallowed = set()
  if exclude is not None:
    for line in exclude:
      disallowed.add( line.strip().upper() )
    log.write( '%i excluded genes added\n' % len(disallowed) )

  # filter the reference
  filtered = 0
  allowed = 0
  found = set()
  blocked = set()
  for line in bed_in:
    if line.startswith('#'):
      bed_out.write( line )
      allowed += 1
    else:
      fields = line.strip().split( '\t' )
      if len(fields) > 3:
        candidate = fields[3].upper()
        if candidate in genes:
          if candidate in disallowed:
            blocked.add( candidate )
          else:
            bed_out.write( line )
            allowed += 1
            found.add( candidate )
        else:
          filtered += 1
      else:
        bed_out.write( line )
        allowed += 1
  
  log.write( '%i lines written, %i lines filtered, %i out of %i candidate genes found\n' % ( allowed, filtered, len(found), len(genes) ) )
  if len(found) != len(genes):
    log.write( 'Not found: %s\n' % ' '.join( sorted( list( genes.difference( found.union( blocked ) ) ) ) ) )
  if len(blocked) > 0:
    log.write( 'Excluded: %s\n' % ' '.join( sorted( list( blocked ) ) ) )

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Convert genelists to bed file')
  parser.add_argument('--exclude', dest='exclude', required=False, help='file containing genes to exclude')
  parser.add_argument('genelists', nargs='*', help='gene list files to include' )
  args = parser.parse_args()
  if args.exclude:
    filter_bed( args.genelists, sys.stdin, sys.stdout, sys.stderr, exclude=open( args.exclude, 'r' ) )
  else:
    filter_bed( args.genelists, sys.stdin, sys.stdout, sys.stderr )

