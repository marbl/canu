
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
/*********************************************************************
   Module:  AS_OVL
   Description:  Assembly Overlap Script Module.  Creates a script
      of overlap jobs that can be run on independent processors
      for large overlap jobs.
   Assumptions:  Input meets specifications in the ProtoIO documents
 *********************************************************************/

/* RCS info
 * $Id: MakeScriptOVL.h,v 1.3 2005-03-22 19:07:07 jason_miller Exp $
 * $Revision: 1.3 $
*/


#ifndef MAKE_SCRIPT_OVL_H
#define MAKE_SCRIPT_OVL_H

//
// Component:
//   make-ovl-script
//
//   Art. Delcher
//   Last revised:  7 Dec 99
//
//   Reference:  "Overlap Module.rtf"   // Not yet updated
// 
//  Description:
// 
//   This component takes a command like an overlapper command
//   and generates a script.  The script is a series of overlap
//   jobs that can be run on separate processors.
// 
// Design:
// 
//   The sizes of the old frag store (old_frag_ct) and the new fragment
//   set (new_frag_ct) are determined.  Then a series of overlap jobs
//   is created.  If the old_frag_ct and new_frag_ct are both sufficiently
//   small, then a single overlap job is emitted.  Otherwise, the
//   first job just commits all the fragments to the store.  Then
//   a series of jobs that overlaps approximately  OLD_SIZE  fragments
//   from the frag store with  NEW_SIZE  fragments in hash tables
//   is generated.  Each job emits a separate .ovl file.  A final
//   step is to concatenate all these .ovl files.
//   
// Limitations:
// 
//   Only the entire old store can be matched against the entire
//   new fragment file.
// 
// Status:
// 
//   Initial implementation.
// 
// Architecture and Dependencies:
// 
//   AS_OVL_delcher.h   Generic  #includes  and simple function prototypes
//   AS_OVL_delcher.c   Corresponding function definitions
//   MakeScriptOVL.h    header file
//   MakeScriptOVL.c    main program
//


#define  MAX_NAME_LEN            500
    //  Longest file name allowed
#define  SCRIPT_NAME             "lsf-ovl"
    //  Default name of script produced by  make-ovl-script
#define  INPUT_FILENAME_EXTENSION   ".urc"
#define  OUTPUT_FILENAME_EXTENSION  ".ovl"
    //  For proto-I/O files

int  OVL_Max_int
    (int, int);
int  OVL_Min_int
    (int, int);
#endif

