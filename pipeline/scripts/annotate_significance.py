#!/usr/bin/env python
# vim: expandtab:ts=4:sw=4:cindent
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

####################################################################################
#
# Purpose:
# 
# This script adds annotations to an Annovar output to rate the
# clinical significance of each variant as a "priority index". These rules are 
# defined by the Bioinformatics working group and are currently 
# specified as follows:
#
#     (1)     Index 1 - missense
#     (1)     Index 2 - rare missense
#             Variant has MAF <0.01 in dbSNP or ESP or 1000G databases
#     (2)     Index 3 - novel missense
#             Variant is not in dbSNP or ESP or 1000G databases 
#     (3)     Index 4 - highly conserved missense
#             Variant has condel score > 0.7 (other possible filters - for discussion)
#     (4)     Index 5 - truncating
#             Out of frame indel (non-recurrent*), nonsense, splice site +/-2bp
#
#     NOTE: these rules were updated 27/5/2014:
#
#             - include 'very rare' variants <0.0005 MAF in cat 2,3
#             - include truncating variants into priority 1 if rare (but not novel)
#
#     NOTE: these rules were updated 2015-6-10:
#
#             - all priorities shifted up by 1
#             - new priority 1 introduced to capture all misense variants 
#                (even non-rare)
#
# Author:   Simon Sadedin, simon.sadedin@mcri.edu.au
# Date:     23/1/2014
#
####################################################################################

import csv
import getopt
import logging as log
import re
import sys

log.basicConfig(level=log.INFO)

class AnnovarPriority:
    """
    Helper class to map Annovar column names to their fields parsed from CSV,
    and to implement logic surrounding categorization.
    """

    # Default MAF threshold for considering a variant 'rare'
    MAF_THRESHOLD = 0.01

    # Default MAF threshold for considering a variant 'very rare'
    MAF_THRESHOLD_VERY_RARE = 0.0005

    # Default Condel Threshold
    CONDEL_THRESHOLD = 0.7

    # Categories of variants as specified by Annovar, mapped to functional categories
    # defined for Melbourne Genomics
    ANNOVAR_EXONIC_FUNCS = {
        "truncating" : set( ( "frameshift insertion","frameshift deletion","frameshift substitution","stopgain SNV","stoploss SNV","stoploss","stopgain", "frameshift_insertion","frameshift_deletion","frameshift_substitution","stopgain_SNV","stoploss_SNV", ) ),
        "missense" : set( ( "nonframeshift insertion","nonframeshift deletion","nonframeshift substitution","nonsynonymous SNV", "nonframeshift_insertion","nonframeshift_deletion","nonframeshift_substitution","nonsynonymous_SNV", ) ),
        "synonymous" : set( ( "synonymous SNV", "synonymous_SNV" ) ),
        "noncoding" : set( ( "intronic","intergenic","ncRNA_intronic","ncRNA_exonic","upstream","downstream","UTR5","UTR3","ncRNA_splicing","upstream;downstream","upstream\x3bdownstream"  ) )
    }

    # These are the Annovar fields that contain population frequency estimates
    # Note we do some fooling around in the maf_value() method to maintain 
    # compatibility with different versions of Annovar
    POPULATION_FREQ_FIELDS = ["esp6500siv2_all", "1000g2014oct_all","exac03"]

    def priority(self):
        """
            Main logic describing how to map any given variant to a clinical significance
            priority index. See the main header for the definition of these categories.

            Note: unknown categories are returned as 9 - that is, extremely high.
        """

        if self.is_missense():
           if self.is_rare():
               if self.is_novel() or self.is_very_rare():
                   if self.is_conserved():
                       return 4 # Missense, novel and conserved => category 4
                   else:
                       return 3 # Missesnse, novel but not highly conserved => category 3
               else:
                    return 2 # Missense & rare but not novel => category 2
           else:
               log.debug("%s:%s is missense and not rare" % (self.Chr,self.Start))
               return 1 # Missense but not even rare => category 1

        elif self.is_truncating():
            # From Natalie, 27/5/2014:
            # With regard to priority 5 truncating variants:
            #  novel should stay in priority 5
            #  rare should be priority 2
            if self.is_novel():
                return 5
            elif self.is_rare():
                log.debug("%s:%s is truncating and rare" % (self.Chr,self.Start))
                return 2
            else:
                log.debug("%s:%s is truncating and not rare" % (self.Chr,self.Start))
                return 1

        elif self.is_noncoding():
            return 0

        elif self.ExonicFunc in ["synonymous SNV", "unknown"]:
            return 0
        else:
            print >>sys.stderr, "WARNING: variant %s:%s %s/%s func=%s failed to be categorized" % \
                    (self.Chr, self.Start, self.Ref, self.Alt, self.ExonicFunc)
            return 9
        
    def is_noncoding(self):
        return self.Func in self.ANNOVAR_EXONIC_FUNCS["noncoding"]

    def is_missense(self):
        log.debug( 'ExonicFunc: %s' % ( self.ExonicFunc ) )
        return self.ExonicFunc in self.ANNOVAR_EXONIC_FUNCS["missense"]

    def is_truncating(self):
        return self.ExonicFunc in self.ANNOVAR_EXONIC_FUNCS["truncating"] or self.Func in ["splicing","exonic;splicing"]

    def is_rare(self):
        # Return true iff at least one database has the variant at > the MAF_THRESHOLD
        log.debug("MAF values for %s:%s are %s", self.Chr, self.Start, map(lambda f: self.maf_value(f),self.POPULATION_FREQ_FIELDS))
        return not any(map(lambda f: self.maf_value(f)>self.MAF_THRESHOLD, self.POPULATION_FREQ_FIELDS))

    def is_very_rare(self):
        # Return true iff at least one database has the variant at > the MAF_THRESHOLD_VERY_RARE
        return not any(map(lambda f: self.maf_value(f)>self.MAF_THRESHOLD_VERY_RARE, self.POPULATION_FREQ_FIELDS))

    def is_novel(self):
        # return true iff the variant has no MAF in any database AND no DBSNP ID
        return not any(map(lambda f: self.maf_value(f) > 0.0, self.POPULATION_FREQ_FIELDS)) and (self.snp138 in ["","."])

    def is_conserved(self):
        # Clarification 27/5/2014:
        # ONLY if condel score is missing, then it can categorised as a 3 if CONSERVED by Annovar
        condel_str = self.Condel 
        if condel_str != "":
            return float(condel_str) >= 0.7
        else:
            return self.phastConsElements46way != ""

class AnnovarLineCSV (AnnovarPriority):
    """
    The function init_columns() must be called, passing the header row as 
    returned by csv.reader() to initialise the class before use.
    """
    # Column names of Annovar file
    columns = []

    def __init__(self, line):
        self.line = line

    @staticmethod
    def init_columns(cols):
        AnnovarLineCSV.columns = cols #+ ["MapQ","QD"]

    def __getattr__(self,name):
        return self.line[self.columns.index(name)]

    def write_row(self, output):
        while len(self.line)<len(AnnovarLineCSV.columns):
            self.line.append("")
        output.writerow(self.line + [av.priority()])
 
    def set_value(self,name,value):
        self.line[self.columns.index(name)]=value

    def maf_value(self, name):
        # Trying to be compatible with multiple versions of Annovar, each having different
        # names for this column
        if name == "exac03" and "ExAC_Freq" in self.columns:
            name = "ExAC_Freq"
        if name == "exac03" and "ExAC_ALL" in self.columns:
            name = "ExAC_ALL"
        value = self.line[self.columns.index(name)]
        if value == "" or value ==".":
            return 0
        else:
            return float(value)


class AnnovarLineVCF (AnnovarPriority):
    FIELDS = { 'Chr': 0, 'Start': 1, 'Id': 2, 'Ref': 3, 'Alt': 4, 'Qual': 5, 'Filter': 6, 'Info': 7 }

    def __init__(self, line):
        self.line = line
        self.fields = line.split( '\t' )
    
    def __getattr__(self,name):
        if name in AnnovarLineVCF.FIELDS:
            return self.fields[AnnovarLineVCF.FIELDS[name]]
        else:
            # key1=val1;key2=vale => { 'key1': 'val1', 'key2': 'val2' }
            info = dict( map( lambda(x): x.split("=",1), [ y for y in self.fields[AnnovarLineVCF.FIELDS['Info']].split( ';' ) if y.find('=') != -1 ] ) )
            if name in info:
                return info[name]
            else:
                return ''

    def maf_value(self, name):
        if name == "exac03" and self.ExAC_Freq != '':
            value = self.ExAC_Freq
        elif name == "exac03" and self.ExAC_ALL != '':
            value = self.ExAC_ALL
        else:
            value = self.__getattr__(name)
        if value == "" or value ==".":
            return 0
        else:
            return float(value)
    
    def write_row( self, output ):
        # add priority to info
        updated_line = re.sub( ';ALLELE_END', ';Priority=%i;ALLELE_END' % self.priority(), self.line )
        output.write( '%s\n' % updated_line )

class AnnovarCSV:
    def __init__( self ):
        self.first = True

    def process( self, line, fh ):
        if self.first:
            self.first = False
            if "Qual" not in line:
                line = line + ["Qual"]
            
            if "Depth" not in line:
                line = line + ["Depth"]
    
            AnnovarLineCSV.init_columns(line)
            header_out = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONE)
            header_out.writerow(AnnovarLineCSV.columns + ["Priority_Index"])
            fh.flush()
            self.output = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        else:
            av = AnnovarLineCSV(line)
            av.write_row(self.output) 

class AnnovarVCF:
    def __init__( self ):
        pass

    def process( self, line, fh ):
        if line.startswith( '#' ):
            fh.write( '%s\n' % line )
        else:
            av = AnnovarLineVCF(line)
            av.write_row( fh ) 

####################################################################################
#
# Main body
#
####################################################################################

def main():
    # Parse command line options
    # annovar file, maf threshold, condel threshold, maf threshold very rare, annovar is vcf (default to csv)
    optstring = "a:f:c:r:v" 
    opts,args = getopt.getopt(sys.argv[1:],optstring)
    
    options = {}
    for opt in opts:
       options[opt[0]] = opt[1]
    
    def usage(msg):
        print >>sys.stderr, "\nERROR: %s\n\nUsage: annotate_significance.py -a <annovar file>\n" % msg
        sys.exit(1)
            
    if not '-a' in options:
        usage("Please provide -a option.")
    
    use_vcf = '-v' in options
    # Note: Annovar does not seem to provide Qual and Depth headings itself
    if '-f' in options:
        AnnovarPriority.MAF_THRESHOLD = float(options['-f'])
    if '-c' in options:
        AnnovarPriority.CONDEL_THRESHOLD = float(options['-c'])
    if '-r' in options:
        AnnovarPriority.MAF_THRESHOLD_VERY_RARE = float(options['-r'])

    # Read the file
    if use_vcf:
        handler = AnnovarVCF()
        reader = open( options['-a'], 'r' ) 
    else:
        handler = AnnovarCSV()
        reader = csv.reader(open(options['-a'], 'r'), delimiter='\t', quotechar='"', doublequote=True)

    # Open CSV writer to standard output, first for header (for body comes in the loop below)
    for line in reader:
        handler.process( line, sys.stdout )
   
if __name__ == "__main__":    
    main()
