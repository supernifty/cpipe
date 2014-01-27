import unittest
import csv

from annotate_significance import Annovar 

class AnnotateSignificanceTest(unittest.TestCase):
    
    
    header = "Func,Gene,ExonicFunc,AAChange,Conserved,SegDup,ESP5400_ALL,1000g2010nov_ALL,dbSNP138,AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,Chr,Start,End,Ref,Obs,Otherinfo,Qual,Depth,Condel".split(",")
    
    # A prototype line that we use to create test data
    line = '"exonic","SCN5A","nonsynonymous SNV","NM_000335:c.G1339T:p.A447S","437;Name=lod=80",,,,,0.42,,,,,,,,,,,,chr3,38646399,38646399,C,A,"het","14.91","19","0.3"\n'
    
    def testNovel(self):
        
        Annovar.init_columns(self.header)
        
        variant = csv.reader([self.line], delimiter=",", quotechar='"').next()
        a = Annovar(variant)
        
        # It has no MAF and no dbSNP ID => is novel
        assert a.is_novel()
        
        # But if we give it a dbSNP ID, then it is not novel
        a.set_value('dbSNP138','rs12345')
        assert not a.is_novel()
        a.set_value('dbSNP138','')
        
        assert a.is_novel()
        
        # Give it  1000 Genomes MAF and it should be not novel
        a.set_value('1000g2010nov_ALL','0.00001')
        assert not a.is_novel()
        
        
    def testPriorityIndex(self):
        
        Annovar.init_columns(self.header)
        
        variant = csv.reader([self.line], delimiter=",", quotechar='"').next()
        a = Annovar(variant)
        
        # Novel, missense, but no Condel score / conserved: category 2
        print a.category()
        assert a.category() == 2
        
        a.set_value('Condel','0.9')
        assert a.category() == 3
        
        a.set_value('Conserved','')
        assert a.category() == 2
        
    
             