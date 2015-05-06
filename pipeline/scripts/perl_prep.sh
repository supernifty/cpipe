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
# This script installs necessary perl modules required by Cpipe (actually, 
# by VEP)
#
###########################################################################

curl -L http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib   
export PERL_MM_OPT="INSTALL_BASE=/home/simon/perl5"
curl -L http://cpanmin.us | perl - Compress::Zlib
curl -L http://cpanmin.us | perl - DBD::mysql
curl -L http://cpanmin.us | perl - DBI
curl -L http://cpanmin.us | perl - Archive::Tar
echo "export PERL5LIB=$HOME/perl5/lib/perl5" >> ~/.bashrc
