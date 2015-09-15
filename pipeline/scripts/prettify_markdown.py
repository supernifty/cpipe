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
#   Turn output from markdown2.py into something that looks better
# Usage:
#   prettify < html > html
####################################################################################

import sys

sys.stdout.write( '<html>\n<head>\n<link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">\n</head>\n<body>\n<div class="pure-g">\n<div class="pure-u-1">\n' )

for line in sys.stdin:
  line = line.replace( '<table>', '<table class="pure-table pure-table-bordered">' )
  sys.stdout.write( line )

sys.stdout.write( '</div>\n</div>\n</body>\n</html>\n' )

