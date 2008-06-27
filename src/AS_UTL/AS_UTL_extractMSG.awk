#=======================================================================
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
#  Awk script for extracting a Celio message of a given type matching
#  a given identifying search target
#
#  Generally used through the shell script extractMSG as follows:
#
#        extractMSG <3-code> <search-target> <input_file>  > output_file
#
#  e.g.  extractMSG IUM acc:4575 < a004.cns > 4575.ium
#
#  Author: Karin A. Remington
#=======================================================================
BEGIN {
        if (length(type) != 3 || search == "" ) {
          print "Usage: "
          print "extractMSG <3-code> <search-target> <input_file> > output_file\n"
          exit
        }
        inmesg = 0;                # flag to indicate whether within an appropriate mesg
        mesgtype = "\{" type "$";  # search string for start of appropriate mesg
        searchline = search "$";   # search string to match candidate messages against
      }
      {
        if (match($0,mesgtype)) {  # start of candidate message
          mesg = $0;
          inmesg = 1;
          ismatch = 0;
          oc = 1;
          next;
        }
        if (inmesg) {              # within a candidate message
          if (match($0,searchline)) {
            ismatch = 1;
            mesg = mesg "\n" $0;
            next;
          }
          if (match($0,"{"))  {
            oc++;                  # increment open braces
            mesg = mesg "\n" $0;   # append current line to message buffer
            next;
          }
          if (match($0,"}"))  {
            oc--;                  # decrement open braces
            mesg = mesg "\n" $0;   # append current line to message buffer
            if (oc == 0) {         # if all open braces have been closed,
              if (ismatch) {       # if match on search string, print
                print mesg;
                # exit;
              }
              inmesg = 0;          # reset
              ismatch = 0;
              mesg = "";
            }
            next;
          }
          mesg = mesg "\n" $0;
        }
     }
END  { }
