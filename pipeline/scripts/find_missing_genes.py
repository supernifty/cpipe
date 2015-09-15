#!/usr/bin/env python2.7
####################################################################################
#
# Melbourne Genomics Batch Validation
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
# * Find genes missing from a bed file
#
# Usage:
#   find_missing_genes bed_file < genes
####################################################################################

import sys

# read bed genes
ref = set()
for line in open( sys.argv[1], 'r' ):
  fields = line.strip().split( '\t' )
  if len(fields) > 3:
    ref.add( fields[3].upper() )

# read genes
missing = set()
for line in sys.stdin:
  if line.startswith( '#' ):
    continue
  fields = line.strip().split('\t')
  candidate = fields[0].upper() 
  if candidate not in ref and candidate != '1' and candidate != 'HGNC_SYMBOL':
    missing.add( candidate )

print '\n'.join( sorted( list( missing ) ) )
