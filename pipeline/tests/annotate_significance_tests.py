# vim: ts=4:expandtab:sw=4:cindent
###########################################################################
#
# This file is part of Cpipe.
# 
# Cpipe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, under version 3 of the License, subject
# to additional terms compatible with the GNU General Public License version 3,
# specified in the LICENSE file that is part of the Cpipe distribution.
#
# Cpipe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################
import unittest
import csv
import re

#
# NOTE: to run these tests use  :
#
#    python annotate_significance_tests.py
#

from annotate_significance import AnnovarLineCSV, AnnovarLineVCF, AnnovarPriority

class AnnotateSignificanceTest(unittest.TestCase):
    
    header = "Func,Gene,ExonicFunc,AAChange,Conserved,SegDup,esp6500siv2_all,1000g2014oct_all,snp138,AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,Chr,Start,End,Ref,Obs,Otherinfo,Qual,Depth,Condel,exac03,phastConsElements46way".split(",")
    
    # A prototype line that we use to create test data
    line_csv = '"exonic","SCN5A","nonsynonymous SNV","NM_000335:c.G1339T:p.A447S","437;Name=lod=80",,,,,0.42,,,,,,,,,,,,chr3,38646399,38646399,C,A,"het","14.91","19","0.3",".",""\n'

    def __init__(self,n):
        unittest.TestCase.__init__(self,n)
        AnnovarLineCSV.init_columns(self.header)
        variant = csv.reader([self.line_csv], delimiter=",", quotechar='"').next()
        self.a = AnnovarLineCSV(variant)
        #print "init ", self.a.columns

    def testNovel(self):
        
        a = self.a
       
        # It has no MAF and no dbSNP ID => is novel
        assert a.is_novel()
        
        # But if we give it a dbSNP ID, then it is not novel
        a.set_value('snp138','rs12345')
        assert not a.is_novel()
        a.set_value('snp138','')
        
        assert a.is_novel()
        
        # Give it  1000 Genomes MAF and it should be not novel
        a.set_value('1000g2014oct_all','0.00001')
        assert not a.is_novel()

        
    def testPriorityIndex(self):
        
        a = self.a
        
        # Novel, missense, no Condel score but conserved by conserved region annotation: category 3
        a.set_value('Condel','')
        a.set_value('phastConsElements46way','blah')
        #print a.priority()
        assert a.priority() == 4
        
        # Adding high Condel score makes it remain priority 3
        a.set_value('Condel','0.9')
        assert a.priority() == 4
        
        # Removing conserved annotation should not make a difference, because there is still a Condel score
        a.set_value('Conserved','')
        assert a.priority() == 4

        # Setting value to < 0.7 should make it not conserved any more, dropping it down to priority 2
        a.set_value('Condel','0.5')
        assert a.priority() == 3

        # Set a conservation value, but it should NOT change because there is a Condel score
        a.set_value('Conserved','437;Name=lod=80')
        assert a.priority() == 3


    def testNovelRare(self):

        # test: variant index 3 can only happen to a novel variant

        a = self.a
        a.set_value('esp6500siv2_all','0.001')
        a.set_value('1000g2014oct_all','0.001')
        a.set_value('Condel','0.9')

        assert a.priority() == 2
        
    def testRareTruncating(self):
        a = self.a
        a.ExonicFunc = "stopgain SNV"
        assert a.priority() == 5

        a.set_value('esp6500siv2_all','0.001')
        assert a.priority() == 2
        
    def testNovelTruncating(self):
        a = self.a
        a.ExonicFunc = "stopgain SNV"

        # Make it novel
        a.set_value('esp6500siv2_all','')
        a.set_value('1000g2014oct_all','')
        a.set_value('snp138','')
        assert a.priority() == 5

    def testMissenseVaryRare(self):
        a = self.a
        a.set_value('esp6500siv2_all','0.0001')
        a.set_value('1000g2014oct_all','')
        a.set_value('Conserved','')
        a.set_value('Condel','')
        assert a.priority() == 3

    def testMissenseCommon(self):
        a = self.a
        a.set_value('esp6500siv2_all','0.5')
        a.set_value('1000g2014oct_all','')
        a.set_value('Conserved','')
        a.set_value('Condel','')
        assert a.priority() == 1

class AnnotateSignificanceVCFTest(unittest.TestCase):


    def __init__(self,n):
        unittest.TestCase.__init__(self,n)
        line_vcf = 'chr3\t38646399\trs2465128\tA\tG\t938.77\t.\tABHet=0.526;AC=1;AF=0.500;AN=2;BaseQRankSum=1.078;DB;DP=98;Dels=0.00;FS=6.746;HaplotypeScore=5.3526;MLEAC=1;MLEAF=0.500;MQ=60.00;MQ0=0;MQRankSum=-0.647;OND=0.010;QD=9.58;ReadPosRankSum=-1.611;CSQ=G|ENSG00000188157|ENST00000379370|Transcript|synonymous_variant|3116|3066|1022|S|tcA/tcG|rs2465128&COSM1126908|ENST00000379370.2:c.3066A>G|ENST00000379370.2:c.3066A>G(p.%3D)||0.45|0.1|0.04|0.09||||ENSP00000368678|0.473157|0.0901792||AGRN|HGNC|YES|;ANNOVAR_DATE=2015-03-22;Func=exonic;Gene=AGRN;GeneDetail=.;ExonicFunc=nonsynonymous_SNV;AAChange=AGRN:NM_198576:exon18:c.A3066G:p.S1022S;genomicSuperDups=.;ExAC_ALL=.;ExAC_AFR=0.5130;ExAC_AMR=0.9312;ExAC_EAS=0.9649;ExAC_FIN=0.9266;ExAC_NFE=0.9018;ExAC_OTH=0.8983;ExAC_SAS=0.8890;snp138=rs2465128;SIFT_score=.;SIFT_pred=.;Polyphen2_HDIV_score=.;Polyphen2_HDIV_pred=.;Polyphen2_HVAR_score=.;Polyphen2_HVAR_pred=.;LRT_score=.;LRT_pred=.;MutationTaster_score=.;MutationTaster_pred=.;MutationAssessor_score=.;MutationAssessor_pred=.;FATHMM_score=.;FATHMM_pred=.;RadialSVM_score=.;RadialSVM_pred=.;LR_score=.;LR_pred=.;VEST3_score=.;CADD_raw=.;CADD_phred=.;GERP++_RS=.;phyloP46way_placental=.;phyloP100way_vertebrate=.;SiPhy_29way_logOdds=.;ALLELE_END GT:AD:DP:GQ:PL  0/1:51,46:98:99:967,0,1036'
        self.a = AnnovarLineVCF(line_vcf)

    def set_value( self, a, key, value ):
        if a.line.find( key ) == -1:
            updated = a.line + ';%s=%s' % ( key, value )
        else:
            updated = re.sub( '%s=[^;]*' % key, '%s=%s' % ( key, value ), a.line )
        return AnnovarLineVCF(updated)

    def testNovel(self):
        
        a = self.a
        a = self.set_value(a, 'snp138','')
       
        # It has no MAF and no dbSNP ID => is novel
        assert a.is_novel()
        
        # But if we give it a dbSNP ID, then it is not novel
        a = self.set_value(a, 'snp138','rs12345')
        assert not a.is_novel()
        a = self.set_value(a, 'snp138','')
        
        assert a.is_novel()
        
        # Give it  1000 Genomes MAF and it should be not novel
        a = self.set_value(a, '1000g2014oct_all','0.00001')
        assert not a.is_novel()

        
    def testPriorityIndex(self):
        
        a = self.a
        
        # Novel, missense, no Condel score but conserved by conserved region annotation: category 3
        a = self.set_value(a, 'Condel','')
        a = self.set_value(a, 'phastConsElements46way','blah')
        #print a.priority()
        assert a.priority() == 4
        
        # Adding high Condel score makes it remain priority 3
        a = self.set_value(a, 'Condel','0.9')
        assert a.priority() == 4
        
        # Removing conserved annotation should not make a difference, because there is still a Condel score
        a = self.set_value(a, 'Conserved','')
        assert a.priority() == 4

        # Setting value to < 0.7 should make it not conserved any more, dropping it down to priority 2
        a = self.set_value(a, 'Condel','0.5')
        assert a.priority() == 3

        # Set a conservation value, but it should NOT change because there is a Condel score
        a = self.set_value(a, 'Conserved','437\x3bName=lod=80')
        assert a.priority() == 3


    def testNovelRare(self):

        # test: variant index 3 can only happen to a novel variant

        a = self.a
        a = self.set_value(a, 'esp6500siv2_all','0.001')
        a = self.set_value(a, '1000g2014oct_all','0.001')
        a = self.set_value(a, 'Condel','0.9')

        assert a.priority() == 2
        
    def testRareTruncating(self):
        a = self.a
        a = self.set_value(a, 'ExonicFunc', "stopgain_SNV" )
        a = self.set_value(a, 'snp138','')
        assert a.priority() == 5

        a = self.set_value(a, 'esp6500siv2_all','0.001')
        assert a.priority() == 2
        
    def testNovelTruncating(self):
        a = self.a
        a = self.set_value( a, 'ExonicFunc', "stopgain_SNV")

        # Make it novel
        a = self.set_value(a, 'esp6500siv2_all','')
        a = self.set_value(a, '1000g2014oct_all','')
        a = self.set_value(a, 'snp138','')
        assert a.priority() == 5

    def testMissenseVaryRare(self):
        a = self.a
        a = self.set_value(a, 'esp6500siv2_all','0.0001')
        a = self.set_value(a, '1000g2014oct_all','')
        a = self.set_value(a, 'Conserved','')
        a = self.set_value(a, 'Condel','')
        assert a.priority() == 3

    def testMissenseCommon(self):
        a = self.a
        a = self.set_value(a, 'esp6500siv2_all','0.5')
        a = self.set_value(a, '1000g2014oct_all','')
        a = self.set_value(a, 'Conserved','')
        a = self.set_value(a, 'Condel','')
        assert a.priority() == 1
 
if __name__ == '__main__':
    unittest.main()             
     
 
if __name__ == '__main__':
    unittest.main()             
