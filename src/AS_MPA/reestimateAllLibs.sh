#!/usr/bin/env bash
#
###########################################################################
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received (LICENSE.txt) a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
# $Id: reestimateAllLibs.sh,v 1.4 2008-06-27 06:29:17 brianwalenz Exp $
#

prefix=${1}

if [ -z ${prefix} ] || [ ${prefix} == "bell-style" ] ; then
  echo "Please specify the prefix of the intra & lib files"
  return
fi

for lib in `cut -f 1 -d ' ' ${prefix}Libs.txt`; do
  reestimateLibs -f ${lib}Lengths.txt \
                            -l ${prefix}Libs.txt \
                            -n ${lib}
done
