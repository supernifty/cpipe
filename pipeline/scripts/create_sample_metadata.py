#!/usr/bin/env python
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
# * add a pipeline run ID to the sample metadata file, read from a pipeline run ID file
# * takes incoming meta from stdin and writes to stdout
####################################################################################

import argparse
import os.path
import random
import sys

def generate_id( f ):
  '''
    given a file, reads the current ID, appends to it, and writes it back to the same file.
    if the file doesn't exist, a random ID is generated.
    format of the ID is site_000000000
  '''
  # NOTE! we don't do any file locking. 
  # parallel pipelines could potentially attempt to update the ID simulatenously, resulting in a non-unique ID
  if os.path.isfile( f ):
    fh = open( f, 'r' )
    current = fh.readline().strip()
    fh.close()
    if "_" in current:
      site, run = current.rsplit("_", 1)
      run = int(run) + 1
      new_id = '%s_%09i' % (site, run)
    elif len(current) == 0: # empty file
      run = 0
      new_id = 'site%i_%09i' % (random.randint(0, 1e6), run )
    else: # no _
      new_id = '%s_%09i' % (current, 0)
  else: # no file
    run = 0
    new_id = 'site%i_%09i' % (random.randint(0, 1e6), run )

  fh = open( f, 'w' )
  fh.write( new_id )
  return new_id

def write( src, target, new_id ):
  '''
    reads lines from src and writes to target, appending the pipeline ID at the start as a new column of a tab separated file
  '''
  first = True
  for line in src:
    if first:
      target.write( 'Pipeline_Run_ID\t%s' % line )
      first = False
    else:
      target.write( '%s\t%s' % ( new_id, line ) )

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generate sample metadata file with pipeline ID')
  parser.add_argument('--id', required=True, help='ID file to read/write')
  args = parser.parse_args() 
  new_id = generate_id( args.id )
  write( sys.stdin, sys.stdout, new_id )
