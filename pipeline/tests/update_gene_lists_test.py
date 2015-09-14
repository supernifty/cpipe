
import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import update_gene_lists

class UpdateGenesListsTest(unittest.TestCase):

    def test_find(self):
      with open( 'CS.genes.txt', 'w' ) as current:
        current.write( '#version\nabc\t1\ndef\t2\n' )
      with open( 'CS.add.genes.txt', 'w' ) as current:
        current.write( 'def\nghi\n' )
      log = StringIO.StringIO()
      update_gene_lists.update_gene_lists( '.', '.', log )
      lines = open( 'CS.genes.txt', 'r' ).readlines()
      assert lines[1:] == [ '#notes 1 gene(s) added: GHI\n', 'ABC\t1\n', 'DEF\t2\n', 'GHI\t1\n' ]
