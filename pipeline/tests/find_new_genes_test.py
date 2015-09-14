import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import find_new_genes

class FindNewGenesTest(unittest.TestCase):

    def test_find(self):
      reference_genes = StringIO.StringIO( 'c1\t1\t2\tabc\nc1\t3\t4\tdef\nc1\t5\t6\tghi' )
      excluded_genes = StringIO.StringIO( 'ghi' )
      sample_lines = ['Cohort\tPrioritised_Genes', 'CS\t4:def,ghi,jkl']
      log = StringIO.StringIO()
      result = find_new_genes.generate_new_genes( sample_lines, log, reference_genes, excluded_genes )
      assert list( result['CS']['notfound'] ) == [ 'JKL' ]
      assert list( result['CS']['add'] ) == [ 'DEF' ]
      assert list( result['CS']['addonce'] ) == [ 'GHI' ]
