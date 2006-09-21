
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
/* 	$Id: Globals_CGW.h,v 1.7 2006-09-21 21:34:00 brianwalenz Exp $	 */
#ifndef GLOBALS_CGW_H
#define GLOBALS_CGW_H

#include "AS_CGW_dataTypes.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_timer.h"
#include "AS_ALN_forcns.h"

/* These are the global data structures for the CGW Module */



typedef struct Global_CGW_tag {
  int verbose;
  MesgReader reader;
  MesgWriter writer;
  MesgWriter errorWriter;
  float transQualityCutoff; // quality cutoff for TransChunkEdges
  uint64 maxSequencedbSize; // maximum size of a sequencedb between checkpoints
  uint64 maxSequencedbCacheSize; // maximum size of cache before flushing
  int32 dumpScaffoldSnapshots;
  int32 maxDegree; // maximum edges to keep for 'nonUnique' nodes
  int32 maxDegreeUnique; // maximum edges to keep for 'Unique' nodes
  int outputCalculatedOffsets;
  int saveCheckPoints;
  int repeatRezLevel;
  int  write_rock_log;  // create log files for rocks iff true
  int walkLevel;   
  /* a variable that contains different alternatives of 
     gap walking */
  int stoneLevel;
  /* a variable that contains different alternatives of 
     stone throwing */
  int  write_stone_log;  // create log files for stones iff true
  int walkScaffoldsBiggestFirst;
  int ignoreChaffUnitigs;
  int performCleanupScaffolds;
  int debugLevel;
  int failOn_NoOverlapFound;
  int cgbUniqueCutoff;
  int cgbDefinitelyUniqueCutoff;
  int cgbApplyMicrohetCutoff;
  int  starting_stone_scaffold;
  float cgbMicrohetProb;
  int  annotateUnitigs;
  int  doInterleavedScaffoldMerging;
  FILE *outfp;  // .cgw    frags, unitigs, singleton contigs and degenerate scaffolds
  FILE *outfp1;  // .cgw_contigs   input for post-cgw consensus
  FILE *outfp2;  // .cgw_scaffolds non-degenerate scaffolds
  FILE *timefp;  // .timing
  FILE *stderrc;  // current - initially set to stderr
#ifdef NEVER
  HISTOGRAM *scaffold_unique;   
  HISTOGRAM *scaffold_repeat;
#endif
  Overlap *(*aligner)(char *, char *, int, int, int, double, double, int, CompareOptions); 
  TimerT RecomputeOffsetsTimer;
  TimerT MergeScaffoldsTimer;
  TimerT BuildSEdgesTimer;
  TimerT InputTimer;
  TimerT OutputTimer;
  TimerT OverlapTimer;
  TimerT ChooseChunksTimer;
  TimerT ConsistencyCheckTimer;
  TimerT WalkUpdateTimer;
  TimerT UpdateTimer;
  TimerT GapFillTimer;
  TimerT GapWalkerTimer;
  TimerT StoneThrowingTimer;
  TimerT BccTimer;
  TimerT ConsensusTimer;
  char TempFileName[1024];
  char Input_File_Name[1024];
  char File_Name_Prefix[1024];
  char Output_File_Name[1024];
  char Frag_Store_Name[1024];
  char Gatekeeper_Store_Name[1024];	  
  char OVL_Store_Name[1024];
}Global_CGW;


Global_CGW *GlobalData;



extern Global_CGW *CreateGlobal_CGW(void);
extern void DeleteGlobal_CGW(Global_CGW *);

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
      exit (EXIT_FAILURE);
    }

  return  fp;
}

#ifdef NEVER
void ResetHistograms_CGW(struct Global_CGW_tag *);
#endif
#endif
