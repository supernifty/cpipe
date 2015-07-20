#!/bin/bash
# vim: ts=4:expandtab:sw=4
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
#
#     Runs all unit tests 
#
########################################################

# python unit tests
pushd pipeline/scripts
python -m unittest discover -s ../tests -p '*_tests.py'
popd

# groovy unit tests
eval `sed 's/\/\/.*$//' pipeline/config.groovy`
pushd pipeline/scripts
JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar ../tests/*.groovy
popd
