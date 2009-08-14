
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

#ifndef GLOBALS_CGW_H
#define GLOBALS_CGW_H

static const char *rcsid_GLOBALS_CGW_H = "$Id: Globals_CGW.h,v 1.28 2009-08-14 13:37:06 skoren Exp $";

#include "AS_CGW_dataTypes.h"
#include "AS_MSG_pmesg.h"
#include "AS_ALN_aligners.h"
#include "AS_UTL_Hash.h"

//  These are the global data structures for the CGW Module

typedef struct Global_CGW_tag {
  int verbose;

  uint64 maxSequencedbSize; // maximum size of a sequencedb between checkpoints
  int repeatRezLevel;
  int stoneLevel;  // a variable that contains different alternatives of stone throwing //
  int ignoreChaffUnitigs;
  int performCleanupScaffolds;
  int debugLevel;
  int cgbUniqueCutoff;
  int cgbDefinitelyUniqueCutoff;
  int cgbApplyMicrohetCutoff;
  float cgbMicrohetProb;
  int  doInterleavedScaffoldMerging;
  int  allowDemoteMarkedUnitigs;

  FILE *cgwfp;    // .cgw            frags, unitigs
  FILE *ctgfp;    // .cgw_contigs    all contigs (input for post-cgw consensus)
  FILE *scffp;    // .cgw_scaffolds  all scaffolds
  FILE *stderrc;  // current - initially set to stderr

  char File_Name_Prefix[FILENAME_MAX];
  char Output_File_Name[FILENAME_MAX];
  char Gatekeeper_Store_Name[FILENAME_MAX];
  char OVL_Store_Name[FILENAME_MAX];

  char unitigOverlaps[FILENAME_MAX];
  int closurePlacement;
}Global_CGW;


extern Global_CGW *GlobalData;

extern Global_CGW *CreateGlobal_CGW(void);
extern void        DeleteGlobal_CGW(Global_CGW *);

extern int         SetFileNamePrefix_CGW(Global_CGW *data, char *name);

/****************************************************************************/
static FILE *  File_Open
(const char * Filename, const char * Mode, int exitOnFailure)

     /* Open  Filename  in  Mode  and return a pointer to its control
      *  block.  If fail, print a message and exit. */

{
  FILE  *  fp;

  fp = fopen (Filename, Mode);
  if  (fp == NULL && exitOnFailure)
    {
      fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
      exit (1);
    }

  return  fp;
}

#ifdef NEVER
void ResetHistograms_CGW(struct Global_CGW_tag *);
#endif
#endif
