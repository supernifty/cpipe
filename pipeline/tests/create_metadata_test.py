import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import create_sample_metadata

class CreateMetadataTest(unittest.TestCase):

    def test_generate_id(self):
       fn = 'tmp%i' % random.randint(0, 1e9)
       with open( fn, 'w' ) as fh:
         fh.write( 'site' )
       
       new_id = create_sample_metadata.generate_id( fn )
       os.remove( fn ) 
       assert new_id == 'site_000000000'

    def test_update_samples(self):
       src = StringIO.StringIO( 'h1\th2\nl1\tl2' )
       target = StringIO.StringIO()
       new_id = create_sample_metadata.write( src, target, 'site_123456789' )
       assert target.getvalue() == 'Pipeline_Run_ID\th1\th2\nsite_123456789\tl1\tl2'
      
