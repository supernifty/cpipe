import unittest
import os
import random
import re

class SelectVariantsTest(unittest.TestCase):

    def test_simple(self):
       # merge command
       os.system( "java -Xmx2g -jar ../../tools/gatk/2.8-1-g932cd3a/GenomeAnalysisTK.jar  -R ../../hg19/ucsc.hg19.fasta -T SelectVariants  --variant sample.vcf  -L sample.bed  --interval_padding 10 -o sample.out" )
       data = open( 'sample.out' ).read()
       assert data.find( '3916' ) == -1
       assert data.find( '955597' ) >= 0 # snp
       assert data.find( '103471456' ) >= 0 # del
       assert data.find( '153233991' ) >= 0 # ins

    def test_indel(self):
       os.system( "java -Xmx2g -jar ../../tools/gatk/2.8-1-g932cd3a/GenomeAnalysisTK.jar  -R ../../hg19/ucsc.hg19.fasta -T SelectVariants  --variant sample.vcf  -L sample.bed  --selectTypeToInclude INDEL --interval_padding 10 -o sample-indel.out" )
       data = open( 'sample-indel.out' ).read()
       assert data.find( '3916' ) == -1
       assert data.find( '955597' ) == -1 # snp
       assert data.find( '103471456' ) >= 0 # del
       assert data.find( '153233991' ) >= 0 # ins

    def test_snv(self):
       os.system( "java -Xmx2g -jar ../../tools/gatk/2.8-1-g932cd3a/GenomeAnalysisTK.jar  -R ../../hg19/ucsc.hg19.fasta -T SelectVariants  --variant sample.vcf  -L sample.bed  --selectTypeToInclude SNP --selectTypeToInclude MIXED --selectTypeToInclude MNP --selectTypeToInclude SYMBOLIC --selectTypeToInclude NO_VARIATION --interval_padding 10 -o sample-snv.out" )
       data = open( 'sample-snv.out' ).read()
       assert data.find( '3916' ) == -1
       assert data.find( '955597' ) >= 0 # snp
       assert data.find( '103471456' ) == -1 # del
       assert data.find( '153233991' ) == -1 # ins


if __name__ == '__main__':
    unittest.main()
