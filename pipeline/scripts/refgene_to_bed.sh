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
#   ./refgene_to_bed.sh < refgene.bed > exons.bed
#
####################################################################################
python refgene_to_bed.py | sort -k1,1 -k2,2n > tmp$$
#$BEDTOOLS/bin/bedtools merge -i tmp$$ -nms | python refgene_to_bed.py post
bedtools merge -i tmp$$ -nms | python refgene_to_bed.py post
rm tmp$$
