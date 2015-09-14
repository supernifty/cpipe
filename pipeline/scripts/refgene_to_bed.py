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
#   Converts hg19_refGene.txt to a bed file
#   Use in conjunction with refgene_to_bed.sh
#
####################################################################################
import datetime
import sys

if len(sys.argv) > 1 and sys.argv[1] == 'post': # split out genes
  sys.stdout.write( '#version %s\n' % datetime.datetime.now().strftime("%Y%m%d") )
  for line in sys.stdin:
    if line.startswith( '#' ):
      sys.stdout.write( line )
    else:
      fields = line.strip().split('\t')
      genes = set( fields[3].split(';') )
      for gene in genes:
        sys.stdout.write( '%s\t%s\t%s\t%s\n' % ( fields[0], fields[1], fields[2], gene ) )
else: # convert from refGene
  sys.stdout.write( '#version %s\n' % datetime.datetime.now().strftime("%Y%m%d") )
  for line in sys.stdin:
    fields = line.strip().split( '\t' )
    exons = ( fields[9].split(','), fields[10].split(',') )
    for exon_start, exon_end in zip( *exons ):
      if exon_start != '' and exon_end != '':
        sys.stdout.write( '%s\t%s\t%s\t%s\n' % ( fields[2], exon_start, exon_end, fields[12] ) )
