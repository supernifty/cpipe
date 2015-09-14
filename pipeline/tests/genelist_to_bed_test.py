import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import genelist_to_bed

class GeneListToBedTest(unittest.TestCase):

    def test_exclude(self):
      bed_in = StringIO.StringIO( 'c1\t1\t2\tabc\nc1\t3\t4\tdef\nc1\t5\t6\tghi' )
      bed_out = StringIO.StringIO()
      log = StringIO.StringIO()
      genelist_to_bed.filter_bed( ['sample_genelist.txt'], bed_in, bed_out, log, exclude=['ghi'] )
      assert bed_out.getvalue() == 'c1\t1\t2\tabc\nc1\t3\t4\tdef\n'
      
    def test_notfound(self):
      bed_in = StringIO.StringIO( 'c1\t1\t2\tabc\nc1\t5\t6\tghi' )
      bed_out = StringIO.StringIO()
      log = StringIO.StringIO()
      genelist_to_bed.filter_bed( ['sample_genelist.txt'], bed_in, bed_out, log, exclude=['ghi'] )
      assert bed_out.getvalue() == 'c1\t1\t2\tabc\n'
      
