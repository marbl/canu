#* $Id: parse_raw_frag_stats.awk,v 1.1.1.1 2004-04-14 13:51:20 catmandew Exp $
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

BEGIN { total_inserts = 0;
        total_deletes = 0;
        total_subs = 0;
        total_basepairs = 0;
        #associate fields with contents
        iid = 1;
        uid = 2;
        type = 3;
        clear_bgn = 4;
        clear_end = 5;
        num_non_errors = 5;
        error_start = num_non_errors + 1;
        if ( placed == "" ) {
          print "Need to define \"placed\" in the command line."
          exit;
        }
        #define files:
        open_files(placed);
      }
      {  
        print celera_inserts;
        num_errors  = (NF - num_non_errors)/2;
        num_inserts = 0;
        num_deletes = 0;
        num_subs = 0;
        frag_len = $(clear_end) - $(clear_bgn);
        if ( $type == "R" ) { 
           insert_file = celera_insert;
           delete_file = celera_delete;
           sub_file = celera_sub;
           all_file = celera_all;
        } else {
           insert_file = ext_insert;
           delete_file = ext_delete;
           sub_file = ext_sub;
           all_file = ext_all;
        }
        for (i=0; i<num_errors; i=i+1) {
           apos = $(2*i+ error_start);
           error_type = $(2*i+ error_start +1);
           rpos = int(apos/frag_len*100);
           if ( error_type == "I" ) {
             num_inserts = num_inserts + 1;
             print rpos >> insert_file
           } else if ( error_type == "D" ) {
             num_deletes = num_deletes + 1;
             print rpos >> delete_file
           } else if ( error_type == "S" ) {
             num_subs = num_subs + 1;
             print rpos >> sub_file
           }
           print rpos >> all_file
        }
        total_inserts += num_inserts;
        total_deletes += num_deletes;
        total_subs += num_subs;
        total_basepairs += frag_len;
      }
END {
    if ( placed == "" ) {
       exit;
    }
    total_errors = total_inserts + total_deletes + total_subs;
    print " Total errors: " total_errors " in " total_basepairs " bp ( " total_errors/total_basepairs*100 " % )";
    print "---------------------------- ";
    print "      inserts: " total_inserts;
    print "      deletes: " total_deletes;
    print "substitutions: " total_subs;
    }

function open_files( placed ) {
    if (placed == "y" ) {
      file_postfix = "placed";
    } else {
      file_postfix = "dregs";
    }
    file_prefix = file_postfix "_celera";
    celera_insert = file_prefix "_insert_errors.cgm";
    celagram_title = "rel. pos. of insertion errors in " file_postfix " Celera Read data";
    print celagram_title >> celera_insert
    celera_delete = file_prefix "_delete_errors.cgm";
    celagram_title = "rel. pos. of deletion errors in " file_postfix " Celera Read data";
    print celagram_title >> celera_delete
    celera_sub = file_prefix "_sub_errors.cgm";
    celagram_title = "rel. pos. of subst. errors in " file_postfix " Celera Read data";
    print celagram_title >> celera_sub
    celera_all = file_prefix "_all_errors.cgm";
    celagram_title = "rel. pos. of all errors in " file_postfix " Celera Read data";
    print celagram_title >> celera_all

    file_prefix = file_postfix "_external";
    ext_insert = file_prefix "_insert_errors.cgm";
    celagram_title = "rel. pos. of insertion errors in " file_postfix " Other data";
    print celagram_title >> ext_insert
    ext_delete = file_prefix "_delete_errors.cgm";
    celagram_title = "rel. pos. of deletion errors in " file_postfix " Other data";
    print celagram_title >> ext_delete
    ext_sub = file_prefix "_sub_errors.cgm";
    celagram_title = "rel. pos. of subst. errors in " file_postfix " Other data";
    print celagram_title >> ext_sub
    ext_all = file_prefix "_all_errors.cgm";
    celagram_title = "rel. pos. of all errors in " file_postfix " Other data";
    print celagram_title >> ext_all
}
