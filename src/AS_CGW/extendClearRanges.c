
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
static char CM_ID[] = "$Id: extendClearRanges.c,v 1.5 2005-04-11 15:09:31 eliv Exp $";


/*********************************************************************
 * Module:  AS_CGW_LoadCheckpoint
 * Description:
 *    For use with debugger to query values in a checkpoint
 * 
 *    Reference: 
 *
 *    Command Line Interface:
 *        $ loadcgw checkpointPath
 *
 *       CGBInputFiles: The file with new IUM,OUM, etc records to process. 
 *
 *       Checkpoint File: File named <outputName>.ckp.n
 *
 * 
 *********************************************************************/
//#define DEBUG 1
//#define DEBUG_BUCIS 1
//#define DEBUG_MERGE_SCAF 1
#define FRAG_POSS_SIZE 200000

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"
#include "AS_PER_SafeIO.h"
#include "ChiSquareTest_CGW.h"

#define NUM_STDDEV_CUTOFF 5.0
#define CONTIG_BASES 1000
#define MAX_EXTENDABLE_FRAGS 100

typedef struct extendableFragT
{
  int fragIid;
  int extension;
  int addedBases;
  int basesToNextFrag;
  int fragOnEnd;
  int unitigID;
} extendableFrag;

typedef struct fragPositionsT
{
  int bgn;
  int end;
} fragPositions;

int findFirstFrag( ContigT *contig, int *fragIid,
                   int *extension, int *basesToNextFrag );
int findLastFrag( ContigT *contig, int *fragIid,
                  int *extension, int *basesToNextFrag );
int dumpFragInfo( int fragIid );
int examineGap( ContigT *lcontig, int lFragIid,
                ContigT *rcontig, int rFragIid, 
                int gapNumber, int *ahang,
                int* olapLengthOut, int *bhang, int *currDiffs,
                int *lcontigBasesIntact, int *rcontigBasesIntact,
                int *closedGapDelta,
                int lBasesToNextFrag, int rBasesToNextFrag,
                int *leftFragFlapLength, int *rightFragFlapLength);
void initVarArrays(void);
void SequenceComplement(char *sequence, char *quality);
LengthT FindGapLength( ChunkInstanceT * lchunk,
                       ChunkInstanceT * rchunk, int verbose);
int *my_Unpack_Alignment_AS(OverlapMesg *align);
int alterContigs_old( int numGaps, int *closedGap, int *closedGapDelta,
                      int *lcontigIdGap, int *rcontigIdGap,
                      int *lcontigLength, int *rcontigLength,
                      int *contigValid, int *allContigLengths,
                      int *alteredScaffoldLengths);
int alterContigs( int numGaps, int *closedGap, int *closedGapDelta,
                  int *lcontigIdGap, int *rcontigIdGap,
                  int *lcontigLength, int *rcontigLength,
                  int *closedGapAhang,
                  int *closedGapOlapLength,
                  int *closedGapBhang,
                  int *closedGapLcontigBasesIntact,
                  int *closedGapRcontigBasesIntact,
                  int *contigValid, int *allContigLengths,
                  int *alteredScaffoldLengths);
int findFirstUnitig( ContigT *contig, int *unitigID );
int findLastUnitig( ContigT *contig, int *unitigID );
void printOverlap( Overlap *olap, char *aseq, char *bseq );
void collectContigStats( int *allContigLengths, int *contigValid);
void produceContigStatsConventional(void);
void produceGapStats( int numGaps, int *closedGap,
                      int *closedGapDelta, int *originalGaps);
void produceScaffoldStats( int *alteredScaffoldLengths );
void produceContigStats( int numGaps, int *closedGap, int *closedGapDelta,
                         int *lcontigIdGap, int *rcontigIdGap,
                         int *lcontigLength, int *rcontigLength,
                         int *contigValid, int *allContigLengths,
                         int numDeletedContigs);
int findRightExtendableFrags( ContigT *contig, int *basesToNextFrag,
                              extendableFrag *extFragsArray );
void extendUnitigs( NodeCGW_T *unitig, int fragIid,
                    extendableFrag extFrag, int scaffoldLeftEnd);
void extendContig( ContigT *contig, int extendAEnd);
int findFirstExtendableFrags( ContigT *contig, extendableFrag *extFragsArray );
int findLastExtendableFrags( ContigT *contig, extendableFrag *extFragsArray );
int GetNewUnitigMultiAlign( NodeCGW_T *unitig, fragPositions *fragPoss,
                            int extendedFragIid );
int setCgwClearRange( int fragIid, int frag3pDelta);
// int getAlteredFragPositions( NodeCGW_T *unitig, fragPositions *fragPoss,
//                             int alteredFragIid, int extension );
void getAlteredFragPositions( NodeCGW_T *unitig, fragPositions **fragPoss,
                              int alteredFragIid, int extension );
void bubbleSortIUMs( IntMultiPos *f_list, int numIMPs);
void adjustUnitigCoords( NodeCGW_T *contig );
void DumpContigMultiAlignInfo ( int contigID );
void DumpContigMultiAlignInfoDirect ( MultiAlignT *cma, int contigID );
void DumpContigUngappedOffsets( int contigID );
void DumpContigMultiAlignInfo ( int contigID );
void DumpContigUngappedOffsets( int contigID );
void DumpUnitigMultiAlignInfo ( int unitigID );

void leftShiftIUM( IntMultiPos *f_list, int numFrags, int extendedFragIid);
void rightShiftIUM( IntMultiPos *f_list, int numFrags, int extendedFragIid);
void saveFragAndUnitigData( int lFragIid, int rFragIid );
void restoreFragAndUnitigData( int lFragIid, int rFragIid );
void saveLocalAlignerVariables( void );
void restoreLocalAlignerVariables( void );
void printGapSizes();
void saveDefaultLocalAlignerVariables( void );
int writeEcrCheckpoint( int *numGapsInScaffold,
                        int *numGapsClosedInScaffold,
                        int *numSmallGapsInScaffold,
                        int *numSmallGapsClosedInScaffold,
                        int *numLargeGapsInScaffold,
                        int *numLargeGapsClosedInScaffold);
int loadEcrCheckpoint( int ckptNum, int *numGapsInScaffold,
                       int *numGapsClosedInScaffold,
                       int *numSmallGapsInScaffold,
                       int *numSmallGapsClosedInScaffold,
                       int *numLargeGapsInScaffold,
                       int *numLargeGapsClosedInScaffold);
int extendCgwClearRange( int fragIid, int frag3pDelta);
void SynchUnitigTWithMultiAlignT( NodeCGW_T *unitig );
int revertToCnsClearRange( int fragIid );

int debug = 0;
int totalContigsBaseChange = 0;

int compExtendableFrags( const void *s1, const void *s2)
{
  const extendableFrag * t1 = s1;
  const extendableFrag * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( t1->extension > t2->extension )
	return -1;
  else if ( t1->extension < t2->extension )
	return 1;
  else 
	return 0;
}

int compIMPs( const void *s1, const void *s2)
{
  const IntMultiPos * t1 = s1;
  const IntMultiPos * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( min( t1->position.bgn, t1->position.end ) < min( t2->position.bgn, t2->position.end ) )
	return -1;
  else if ( min( t1->position.bgn, t1->position.end ) > min( t2->position.bgn, t2->position.end ) )
	return 1;
  else 
  {
	if ( t1 > t2 )
	  return -1;
	else if ( t2 < t1 )
	  return 1;
	else  // should never happen
	  return 0;
  }
}

static ReadStructp fsread = NULL;
static VA_TYPE(char) *lContigConsensus = NULL;
static VA_TYPE(char) *rContigConsensus = NULL;
static VA_TYPE(char) *lContigQuality = NULL;
static VA_TYPE(char) *rContigQuality = NULL;
static VA_TYPE(char) *reformed_consensus = NULL;
static VA_TYPE(char) *reformed_quality = NULL;
static VA_TYPE(int32) *reformed_deltas = NULL;
static VA_TYPE(IntElementPos) *ContigPositions = NULL;
static VA_TYPE(IntElementPos) *UnitigPositions = NULL;

// externable variables for controlling use of Local_Overlap_AS_forCNS
extern int MaxGaps;		    // [initialized value is 12 -- no more than this many segments in the chain]
extern int MaxBegGap;       // [ init value is 200; this could be set to the amount you extend the clear 
                            // range of seq b, plus 10 for good measure]
extern int MaxEndGap;	    // [ init value is 200; this could be set to the amount you extend the 
                            // clear range of seq a, plus 10 for good measure]
extern int MaxInteriorGap;	// [ initial value is 1000 (should have almost no effect)
                            // and defines the largest gap between segments in the chain]
                            // Also: allowed size of gap within the alignment
                            // -- forcing relatively good alignments, compared to those
 				            // allowed in bubble-smoothing where indel polymorphisms are expected
extern int asymmetricEnds;  // boolean to cause the size of an "end gap" to
                            // be evaluated with regard to the clear range extension


int main( int argc, char *argv[])
{
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStoreName = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE, singleSid = NULLINDEX;
  int ckptNum = NULLINDEX;
  int i; // , index;
  int sid, startingGap = 0, setStartingGap = FALSE;
  int numExtendableGaps = 0, numGaps = 0, numGapsClosed = 0, totalBasesInClosedGaps = 0;
  int leftContigExtendable = 0, rightContigExtendable = 0, bothContigsExtendable = 0;
  int gapNumber = 0, numSmallGaps = 0, numSmallGapsClosed = 0, numLargeGaps = 0, numLargeGapsClosed = 0;
  int numSmallGapsThisScaff, numSmallGapsClosedThisScaff, numLargeGapsThisScaff, numLargeGapsClosedThisScaff;
  float maxGapSizeClosed = 0.0;
  int maxGapSizeClosedNumber = -1, numGapsVarTooSmall = 0;
  int totalOlapLength = 0, totalOlapDiffs = 0;
  float totalOlapVariance = 0.0;
  int surrogateOnLeftEnd, surrogateOnRightEnd, surrogatesOnEnd = 0;
  time_t t1;
  int *originalGaps, *alteredScaffoldLengths;
  int *closedGap, *closedGapDelta, *lcontigIdGap, *rcontigIdGap, *lcontigLength, *rcontigLength;
  int *contigValid, *allContigLengths;
  // int numDeletedContigs;
  int *closedGapAhang, *closedGapOlapLength, *closedGapBhang;
  int *closedGapLcontigBasesIntact, *closedGapRcontigBasesIntact;
  int numClosingsTried = 0, unitigMultiAlignFailures = 0, replaceEndUnitigFailures = 0, createAContigFailures = 0;
  int noOverlapFound = 0;
  // fragPositions *fragPoss;
  double sumScaffoldLengths = 0, sumScaffoldLengthsLastCkp = 0;
  int *numGapsInScaffold, *numGapsClosedInScaffold;
  int *numSmallGapsInScaffold, *numSmallGapsClosedInScaffold;
  int *numLargeGapsInScaffold, *numLargeGapsClosedInScaffold;
  
  // save off whatever the rest of the world has for default values for Local_Overlap_AS_forCNS
  saveDefaultLocalAlignerVariables();
  
  // set some variables to control Local_Overlap_AS_forCNS
  MaxGaps = 5;
  asymmetricEnds = TRUE;
  MaxInteriorGap = 30;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("extendClearRanges.stderr","w");
  data->timefp = stderr;
  data->logfp = stderr;
  
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:C:f:g:m:n:s:")) != EOF)){
      switch(ch) {
		case 'c':
		{
		  strcpy( data->File_Name_Prefix, argv[optind - 1]);
		  setPrefixName = TRUE;		  
		}
		break;
		case 'C':
		  startingGap = atoi(argv[optind - 1]);
		  setStartingGap = TRUE;
		  break;
		case 'f':
		{
		  strcpy( data->Frag_Store_Name, argv[optind - 1]);
		  setFragStoreName = TRUE;
		}
		break;
		case 'g':
		{
		  strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
		  setGatekeeperStore = TRUE;
		}
		break;	  
		case 'm':
		  MaxInteriorGap = atoi(argv[optind - 1]);
		  fprintf( stderr, "setting MaxInteriorGap to %d\n", MaxInteriorGap);
		  break;
		case 'n':
		  ckptNum = atoi(argv[optind - 1]);
		  break;
		case 's':
		  singleSid = atoi(argv[optind - 1]);
		  setSingleSid = TRUE;
		  fprintf( stderr, "setting singleSid to %d\n", singleSid);
		  break;
		case '?':
		  fprintf(stderr,"Unrecognized option -%c",optopt);
		default :
		  errflg++;
      }
    }
    if((setPrefixName == FALSE) || (setFragStoreName == 0) || (setGatekeeperStore == 0) || (setSingleSid == 0))
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStoreName = %d setGatekeeperStore = %d outputPath = %s "
			"setSingleSid = %d\n",
		argc, optind, setFragStoreName, setGatekeeperStore, outputPath, setSingleSid);
	fprintf (stderr, "USAGE:  loadcgw -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>"
	  "-s <StartingScaffoldNumber>\n");
	exit (EXIT_FAILURE);
      }
  }

  if (setStartingGap == TRUE)
	fprintf( stderr, "set starting gap to %d\n", startingGap);

  t1 = time(0);
  fprintf( stderr, "====> Starting at %s\n", ctime(&t1));

  // check to see if db.frg.orig exists, and if not, copy db.frg to it
  {
	  char temp_buf[1024];
	  int sysReturn;
	  FILE *fileExistenceCheck;
	  
	  sprintf( temp_buf, "%s/db.frg.orig", GlobalData->Frag_Store_Name);
	  fileExistenceCheck = fopen( temp_buf, "r" );
	  if ( fileExistenceCheck == NULL )
	  {
		  fprintf( stderr, "file %s/db.frg.orig does not exist, creating upon start from ckpt %d.\n",
				   GlobalData->Frag_Store_Name, ckptNum );
		  
		  sprintf( temp_buf, "cp %s/db.frg %s/db.frg.orig", 
				   GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
		  sysReturn = system( temp_buf );
		  if ( sysReturn != -1)
			  fprintf( stderr, "copied %s/db.frg to %s/db.frg.orig\n",
					   GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
		  else
		  {
			  fprintf( stderr, "error copying %s/db.frg to %s/db.frg.orig\n",
					   GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
			  assert(0);
		  }
		  
		  sprintf( temp_buf, "chmod 444 %s/db.frg.orig", GlobalData->Frag_Store_Name);
		  sysReturn = system( temp_buf );
		  if ( sysReturn == -1)
		  {
			  fprintf( stderr, "error doing chmod on %s/db.frg.orig\n", GlobalData->Frag_Store_Name);
			  assert(0);
		  }
	  }
  }

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, TRUE);
  GlobalData->aligner=Local_Overlap_AS_forCNS;
  initVarArrays();

  // localeCam();

  closedGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGap == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGap\n");
        fprintf( stderr, "Tried to get " F_SIZE_T " bytes\n",
                 GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	assert(0);
  }
  closedGapDelta = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapDelta == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapDelta\n");
	assert(0);
  }  

  // following arrays are used in checkpointing
  numGapsInScaffold = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  if (numGapsInScaffold == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for numGapsInScaffold\n");
	assert(0);
  }  
  numGapsClosedInScaffold = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  if (numGapsClosedInScaffold == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for numGapsClosedInScaffold\n");
	assert(0);
  }  
  numSmallGapsInScaffold = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  if (numSmallGapsInScaffold == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for numSmallGapsInScaffold\n");
	assert(0);
  }  
  numSmallGapsClosedInScaffold = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  if (numSmallGapsClosedInScaffold == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for numSmallGapsClosedInScaffold\n");
	assert(0);
  }  
  numLargeGapsInScaffold = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  if (numLargeGapsInScaffold == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for numLargeGapsInScaffold\n");
	assert(0);
  }  
  numLargeGapsClosedInScaffold = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  if (numLargeGapsClosedInScaffold == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for numLargeGapsClosedInScaffold\n");
	assert(0);
  }  
  loadEcrCheckpoint( ckptNum, numGapsInScaffold, numGapsClosedInScaffold,
					 numSmallGapsInScaffold, numSmallGapsClosedInScaffold,
					 numLargeGapsInScaffold, numLargeGapsClosedInScaffold);

  for ( i = 0; i < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++)
  {
	numGapsInScaffold[i] = 0;
	numGapsClosedInScaffold[i] = 0;
	numSmallGapsInScaffold[i] = 0;
	numSmallGapsClosedInScaffold[i] = 0;
	numLargeGapsInScaffold[i] = 0;
	numLargeGapsClosedInScaffold[i] = 0;
  }
  

  // we don't use these anymore, they are from when we simulated gap closing
#if 1
  if (1)
  {
	originalGaps = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (originalGaps == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for originalGaps\n");
	  assert(0);
	}
	alteredScaffoldLengths = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
	if (alteredScaffoldLengths == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for alteredScaffoldLengths\n");
	  assert(0);
	}
	lcontigIdGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (lcontigIdGap == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for lcontigIdGap\n");
	  assert(0);
	}
	rcontigIdGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (rcontigIdGap == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for rcontigIdGap\n");
	  assert(0);
	}
	lcontigLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (lcontigLength == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for lcontigLength\n");
	  assert(0);
	}
	rcontigLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (rcontigLength == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for rcontigLength\n");
	  assert(0);
	}
	contigValid = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (contigValid == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for contigValid\n");
	  assert(0);
	}
	allContigLengths = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (allContigLengths == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for allContigLengths\n");
	  assert(0);
	}
	closedGapAhang = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (closedGapAhang == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for closedGapAhang\n");
	  assert(0);
	}
	closedGapOlapLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (closedGapOlapLength == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for closedGapOlapLength\n");
	  assert(0);
	}
	closedGapBhang = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (closedGapBhang == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for closedGapBhang\n");
	  assert(0);
	}
	closedGapLcontigBasesIntact = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (closedGapLcontigBasesIntact == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for closedGapLcontigBasesIntact\n");
	  assert(0);
	}
	closedGapRcontigBasesIntact = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
	if (closedGapRcontigBasesIntact == NULL)
	{
	  fprintf( stderr, "Could not safe_malloc space for closedGapRcontigBasesIntact\n");
	  assert(0);
	}
  }
#endif
  
#if 0
  for ( i = 0; i < GetNumGraphNodes(ScaffoldGraph->ContigGraph); i++)
	contigValid[i] = FALSE;
  collectContigStats( allContigLengths, contigValid);
#endif

  // fprintf( stderr, "before anything happens:\n");
  // DumpContigMultiAlignInfo ( 74812 );

  // print initial separation of contigs
  // fprintf( stderr, "initial gap sizes\n -------------------------------------------------\n");
  // printGapSizes();

  //
  // scan all the scaffolds
  //

  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
        // CIScaffoldTIterator CIsTemp;
	CIScaffoldT * scaff;
	int icnt, lFragIid, rFragIid, lextension, rextension;
	// int leftFragFound, rightFragFound;
	extendableFrag leftExtFragsArray[ MAX_EXTENDABLE_FRAGS ], rightExtFragsArray[ MAX_EXTENDABLE_FRAGS ];
	int numLeftFrags, numRightFrags, rcontigID,lcontigID;
	IntElementPos contigPos;
	ContigT *lcontig, *rcontig, *newContig;
	int lunitigID, runitigID;

	// speed hack
	if ( singleSid != -1 )
	{
	  if ( sid != singleSid )
		continue;
	  else  // we want all the scaffolds from singleSid on
		singleSid = -1;
	}

	
 	lextension = rextension = 0;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    


	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    

	if ((isDeadCIScaffoldT(scaff)) ||
		(scaff->type != REAL_SCAFFOLD) ||
		(scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}


	{
	  time_t tt = time(0);

	  fprintf(stderr,"\n=====================================================================\n");
	  fprintf(stderr,"=== examing scaffold %d, size %f at time %s\n", sid, scaff->bpLength.mean, ctime(&tt));
	}
	numSmallGapsThisScaff = 0;
	numSmallGapsClosedThisScaff = 0;
	numLargeGapsThisScaff = 0;
	numLargeGapsClosedThisScaff = 0;

	// initialize
	// alteredScaffoldLengths[ sid ] = (int) scaff->bpLength.mean;
	// fprintf( stderr, "alteredScaffoldLengths[ %d ] initted to %d\n",
	// 	 sid, alteredScaffoldLengths[ sid ]);


	
	// make sure the scaffold is there
	assert(scaff != NULL);

	icnt = 0;
	// InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
	// while (NextCIScaffoldTIterator(&CIsTemp))

	lcontig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.AEndCI);
    lcontigID = lcontig->id;
	rcontigID = lcontig->BEndNext;
	while ( rcontigID != -1 )
	{
	  NodeOrient lcontigOrientation, rcontigOrientation;
	  LengthT gapSize, newOffsetAEnd, newOffsetBEnd;
	  double maxRContigOffset;
	  int32 nextContigIndex;
	  
#if 0	  
	  // not walking off of scaffolds currently
	  if (CIsTemp.next == -1)
		break;
		
	  //
	  // find the chunks in the gap by walking between the chunk <CIs.curr>
	  // and the chunk <CIs.next>
	  //

	  lcontig = GetGraphNode( ScaffoldGraph->RezGraph, CIsTemp.curr);
	  rcontig = GetGraphNode( ScaffoldGraph->RezGraph, CIsTemp.next);
#endif

#if 1
	  fprintf( stderr, "at top of loop: lcontig->BEndNext: %d\n", lcontig->BEndNext);
	  rcontig = GetGraphNode( ScaffoldGraph->ContigGraph, rcontigID);
#endif

	  // lcontigIdGap[ gapNumber ] = lcontig->id;
	  // rcontigIdGap[ gapNumber ] = rcontig->id;	  
	  // lcontigLength[ gapNumber ] = (int) lcontig->bpLength.mean;
	  // rcontigLength[ gapNumber ] = (int) rcontig->bpLength.mean;

	  assert(lcontig != NULL);
	  assert(rcontig != NULL);		

	  surrogatesOnEnd = 0;
	  gapSize = FindGapLength( lcontig, rcontig, FALSE);
	  // originalGaps[numGaps] = (int) gapSize.mean;
	  numGaps++;

	  if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
		lcontigOrientation = A_B;
	  else
		lcontigOrientation = B_A;

	  if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
		rcontigOrientation = A_B;
	  else
		rcontigOrientation = B_A;

	  fprintf( stderr, "\n\n\n---------------------------------------------------------------\n");
	  {
		time_t tt = time(0);
		fprintf( stderr, "time at gap %8d at %s\n", gapNumber, ctime(&tt));
	  }
	  fprintf( stderr, "examining gap %d from %d (orient: %c, pos: %f, %f) to %d (orient: %c, pos: %f, %f), size: %lf \n", 
			   gapNumber, 
			   lcontig->id, lcontigOrientation, lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean,
			   rcontig->id, rcontigOrientation, rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean,
			   gapSize.mean);
	  
	  if (gapSize.mean < 100.0)
	  {
		numSmallGaps++;
		numSmallGapsThisScaff++;
	  }
	  else
	  {
		numLargeGaps++;
		numLargeGapsThisScaff++;
	  }
		
	  fprintf( stderr, "\nexamining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation);
	  lFragIid = - 1;
	  numLeftFrags = numRightFrags = 0;
	  // find the extreme read on the correct end of the lcontig
	  if ( lcontigOrientation == A_B )
	  {
		numLeftFrags = findLastExtendableFrags( lcontig, leftExtFragsArray );
		// leftFragFound = findLastFrag( lcontig, &lFragIid, &lextension, &basesToNextFrag );
		surrogateOnLeftEnd = findLastUnitig( lcontig, &lunitigID );
	  }
	  else
	  {
		numLeftFrags = findFirstExtendableFrags( lcontig, leftExtFragsArray );
		// leftFragFound = findFirstFrag( lcontig, &lFragIid, &lextension, &basesToNextFrag );
		surrogateOnLeftEnd = findFirstUnitig( lcontig, &lunitigID );
	  }
	  surrogatesOnEnd += surrogateOnLeftEnd;
	  if ( surrogateOnLeftEnd ) numLeftFrags = 0;
#if 0	  
	  if ( numLeftFrags > 0 )
	  {
		leftFragFound = TRUE;
		if (lFragIid != leftExtFragsArray[0].fragIid)
		  assert(0);
		lFragIid = leftExtFragsArray[0].fragIid;
	  }
#endif
	  
	  fprintf( stderr, "\nexamining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation);
	  rFragIid = - 1;
	  // find the extreme read on the correct end of the rchunk
	  if ( rcontigOrientation == A_B )
	  {
		numRightFrags = findFirstExtendableFrags( rcontig, rightExtFragsArray );
		// rightFragFound = findFirstFrag( rcontig, &rFragIid, &rextension, &basesToNextFrag );
		surrogateOnRightEnd = findFirstUnitig( rcontig, &runitigID );
	  }
	  else
	  {
		numRightFrags = findLastExtendableFrags( rcontig, rightExtFragsArray );
		// rightFragFound = findLastFrag( rcontig, &rFragIid, &rextension, &basesToNextFrag );
		surrogateOnRightEnd = findLastUnitig( rcontig, &runitigID );
	  }
	  surrogatesOnEnd += surrogateOnRightEnd;
	  if ( surrogateOnRightEnd ) numRightFrags = 0;
#if 0
	  if ( numRightFrags > 0 )
	  {
		rightFragFound = TRUE;
		if (rFragIid != rightExtFragsArray[0].fragIid)
		  assert(0);
		rFragIid = rightExtFragsArray[0].fragIid;
	  }
#endif

	  if ( setStartingGap == TRUE && gapNumber < startingGap ) numLeftFrags = numRightFrags = 0;

	  // hack 
	  // if ( rcontig->id == 405527 ) numLeftFrags = numRightFrags = 0;


#define DIAG_PRINTS_off
	  
	  if ( numLeftFrags > 0 ) leftContigExtendable++;
	  if ( numRightFrags > 0 ) rightContigExtendable++;
	  if ( numLeftFrags > 0 && numRightFrags > 0 ) bothContigsExtendable++;

	  {	  
		// this test should be applied after the extensions have been factored into gap size
		// original:
		// if ( gapSize.mean - lextension - rextension > NUM_STDDEV_CUTOFF * sqrt(gapSize.variance) && gapSize.mean > 100.0)
		// hacking
		// compute here using {left. right}ExtFragsArray[ index ].extension, not lextension and rextension

		int ahang, currLength, bhang, lcontigBasesIntact, rcontigBasesIntact;
		int currDiffs;
		int leftFragIndex, rightFragIndex;
		int leftFragFlapLength, rightFragFlapLength;
		int gotNewLeftMA, gotNewRightMA;
		InfoByIID *info;
		
		fprintf( stderr, "gap (%d, %d) has extendable frag pair on left: %d, on right: %d\n", 
				 lcontig->id, rcontig->id, numLeftFrags, numRightFrags);
		numExtendableGaps++;
		
		// set an extra member of the arrays for when the contig is not being extended
		leftExtFragsArray[ numLeftFrags ].fragIid = -1;
		leftExtFragsArray[ numLeftFrags++ ].extension = 0;
		rightExtFragsArray[ numRightFrags ].fragIid = -1;		  
		rightExtFragsArray[ numRightFrags++ ].extension = 0;
		closedGap[ gapNumber ] = FALSE;
			
		for (leftFragIndex = 0; leftFragIndex < numLeftFrags && closedGap[ gapNumber ] == FALSE; leftFragIndex++)
		{
		  for (rightFragIndex = 0; rightFragIndex < numRightFrags && closedGap[ gapNumber ] == FALSE; rightFragIndex++)
		  {
            fprintf( stderr, "examining frags %d and %d\n", leftExtFragsArray[ leftFragIndex ].fragIid, 
					 rightExtFragsArray[ rightFragIndex ].fragIid);

			if ( gapSize.mean - leftExtFragsArray[ leftFragIndex ].extension - 
				 rightExtFragsArray[ rightFragIndex ].extension > 
				 NUM_STDDEV_CUTOFF * sqrt(gapSize.variance) && gapSize.mean > 100.0)
			{
			  fprintf( stderr, "leftExtFragsArray[ %d ].extension: %10d, rightExtFragsArray[ %d ].extension: %10d\n",
					   leftFragIndex, leftExtFragsArray[ leftFragIndex ].extension, 
					   rightFragIndex, rightExtFragsArray[ rightFragIndex ].extension);
			  // numGapsVarTooSmall++;
			  closedGap[ gapNumber ] = FALSE;
			  fprintf( stderr, "gap variance too large (gapSize - extensions: %.2f, %.1f * sqrt(gapSize.variance): %.2f\n",
					   gapSize.mean - leftExtFragsArray[ leftFragIndex ].extension - 
					   rightExtFragsArray[ rightFragIndex ].extension, 
					   NUM_STDDEV_CUTOFF, NUM_STDDEV_CUTOFF * sqrt(gapSize.variance));
			  continue;
			}

			lFragIid = leftExtFragsArray[ leftFragIndex ].fragIid;
			rFragIid = rightExtFragsArray[ rightFragIndex ].fragIid;

			// have to check and make sure that the frags belong to the correct unitig
			if ( lFragIid != -1)
			{
			  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, lFragIid);
			  assert( info->set );
			  if (GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex)->cid != lunitigID)
				continue;
			}
			
			if ( rFragIid != -1)
			{
			  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, rFragIid);
			  assert( info->set );
			  if (GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex)->cid != runitigID)
				continue;
			}

			numClosingsTried++;
			  
			// dumpContigInfo( lcontig );
			// dumpContigInfo( rcontig );
			  
			if ( examineGap( lcontig, lFragIid, rcontig, rFragIid, 
							 gapNumber, &ahang, &currLength, &bhang, &currDiffs,
							 &lcontigBasesIntact, &rcontigBasesIntact, &closedGapDelta[gapNumber],
							 leftExtFragsArray[ leftFragIndex ].basesToNextFrag, 
							 rightExtFragsArray[ rightFragIndex ].basesToNextFrag,
							 &leftFragFlapLength, &rightFragFlapLength))
			{
			  int keepGap = TRUE;

			  // save off copies of everything we might alter so we can restore if gap closing fails
			  saveFragAndUnitigData( lFragIid, rFragIid );

			  totalBasesInClosedGaps += (int) gapSize.mean;
				
			  if ( CONTIG_BASES < 2000 )
			  {
				// these checks Granger suggested
				if (ahang + currLength + bhang - 1000 > min( CONTIG_BASES, (int) lcontig->bpLength.mean) +
					min( CONTIG_BASES, (int) rcontig->bpLength.mean))
				{
				  fprintf( stderr, "at gapNumber %d, ahang + currLength + bhang - 1000 = %d\n",
						   gapNumber, ahang + currLength + bhang - 1000);
				  fprintf( stderr, "at gapNumber %d, back(A) + back(B) = %d\n",
						   gapNumber, min( CONTIG_BASES, (int) lcontig->bpLength.mean) +
						   min( CONTIG_BASES, (int) rcontig->bpLength.mean));
				  keepGap = FALSE;
				}
				
				if (ahang + currLength + bhang + 700 < min( CONTIG_BASES, (int) lcontig->bpLength.mean) +
					min( CONTIG_BASES, (int) rcontig->bpLength.mean))
				{
				  fprintf( stderr, "at gapNumber %d, ahang + currLength + bhang + 500 = %d\n",
						   gapNumber, ahang + currLength + bhang + 500);
				  fprintf( stderr, "at gapNumber %d, back(A) + back(B) = %d\n",
						   gapNumber, min( CONTIG_BASES, (int) lcontig->bpLength.mean) +
						   min( CONTIG_BASES, (int) rcontig->bpLength.mean));
				  keepGap = FALSE;
				}
			  }

			  gotNewLeftMA = gotNewRightMA = TRUE;

			  if (keepGap)
			  {
				// extend the clear ranges of the frags
				if (lFragIid != -1)
				{
				  if (leftExtFragsArray[ leftFragIndex ].addedBases - leftFragFlapLength < 0)
					// || leftExtFragsArray[ leftFragIndex ].extension - leftFragFlapLength < 0)
					keepGap = FALSE;
				  else {
				    fprintf(stderr,"adjusting left frg clear range by %d - %d = %d bases\n",
					    leftExtFragsArray[ leftFragIndex ].addedBases,leftFragFlapLength,
					    leftExtFragsArray[ leftFragIndex ].addedBases - leftFragFlapLength);

					extendCgwClearRange( lFragIid,
										 leftExtFragsArray[ leftFragIndex ].addedBases - leftFragFlapLength);
				  }
				}
				if (rFragIid != -1)
				{
				  if (rightExtFragsArray[ rightFragIndex ].addedBases - rightFragFlapLength < 0)
					keepGap = FALSE;
				  else{
				    fprintf(stderr,"adjusting right frg clear range by %d - %d = %d bases\n",
					    rightExtFragsArray[ rightFragIndex ].addedBases,rightFragFlapLength,
					    rightExtFragsArray[ rightFragIndex ].addedBases - rightFragFlapLength);

					extendCgwClearRange( rFragIid,
										 rightExtFragsArray[ rightFragIndex ].addedBases - rightFragFlapLength); 
				  }
				}
			  }

			  if (keepGap)  // the fragment extensions have succeeded
			  {
				// InfoByIID *info;
				CIFragT *frag;
				NodeCGW_T *unitig;
				MultiAlignT *new_cma;
				fragPositions *fragPoss;
				  
				fprintf( stderr, "before anything:\n");
				DumpContigMultiAlignInfo ( lcontig->id );
				DumpContigUngappedOffsets( lcontig->id );
				DumpContigMultiAlignInfo ( rcontig->id );
				DumpContigUngappedOffsets( rcontig->id );

				if (1)
				{
				  if (lcontigOrientation == A_B)
					fprintf( stderr, "before altering, lctg: %12.0f, %12.0f\n",
							 lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean);
				  else
					fprintf( stderr, "before altering, lctg: %12.0f, %12.0f\n",
							 lcontig->offsetBEnd.mean, lcontig->offsetAEnd.mean);
				  if (rcontigOrientation == A_B)
					fprintf( stderr, "before altering, rctg: %12.0f, %12.0f\n",
							 rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean);
				  else
					fprintf( stderr, "before altering, rctg: %12.0f, %12.0f\n",
							 rcontig->offsetBEnd.mean, rcontig->offsetAEnd.mean);
				}

				// save the max offset of the right contig so we know how to adjust the offsets of
				// the contigs further along the scaffold later
				maxRContigOffset = max( rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean );
				nextContigIndex = rcontig->BEndNext;
				  
				// alter the unitigs in cgw memory struct land
					
				// left unitig
				if (lFragIid != -1)
				{
				  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, lFragIid);
				  assert(info->set);
				  frag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
					
				  unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->CIid);
				  // extendUnitigs( unitig, lFragIid, leftExtFragsArray[ leftFragIndex ], TRUE);
				  getAlteredFragPositions( unitig, &fragPoss, lFragIid, 
										   leftExtFragsArray[ leftFragIndex ].extension - leftFragFlapLength );
				  saveLocalAlignerVariables();
				  gotNewLeftMA = GetNewUnitigMultiAlign( unitig, fragPoss, lFragIid );
				  restoreLocalAlignerVariables();
				  free( fragPoss );  // inefficient, let's just do a big array once that get's reused

				  if ( !gotNewLeftMA )
					unitigMultiAlignFailures++;

				  if ( gotNewLeftMA )
				  {
					int extendToLeft;
					SynchUnitigTWithMultiAlignT( unitig );
					
					fprintf( stderr, "before ReplaceEndUnitigInContig:\n");
					DumpContigMultiAlignInfo ( lcontig->id );

					if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
					  extendToLeft = FALSE;
					else
					  extendToLeft = TRUE;

					saveLocalAlignerVariables();
					new_cma = ReplaceEndUnitigInContig( ScaffoldGraph->sequenceDB,
														ScaffoldGraph->fragStore,
														lcontig->id, unitig->id, extendToLeft, GlobalData->aligner);
					restoreLocalAlignerVariables();

					if ( new_cma )
					{
					  UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, lcontig->id, FALSE);
					  InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, lcontig->id, 
													 FALSE, new_cma, FALSE);
					}
					else
					  gotNewLeftMA = FALSE;
					  
					if ( !new_cma )
					  replaceEndUnitigFailures++;

					fprintf( stderr, "strlen( Getchar (new_cma->consensus)): " F_SIZE_T "\n",
							 strlen( Getchar (new_cma->consensus, 0)));
					  
					// fprintf( stderr, "after ReplaceEndUnitigInContig:\n");
					// DumpContigMultiAlignInfo ( lcontig->id );
					  
					if ( lcontigOrientation == A_B )
					  extendContig( lcontig, FALSE);
					else
					  extendContig( lcontig, TRUE);
				  }
				}
				  
				// right unitig
				if (rFragIid != -1)
				{
				  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, rFragIid);
				  assert(info->set);
				  frag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
					
				  unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->CIid);
				  // extendUnitigs( unitig, rFragIid, rightExtFragsArray[ rightFragIndex ], FALSE);
					
				  getAlteredFragPositions( unitig, &fragPoss, rFragIid, 
										   rightExtFragsArray[ rightFragIndex ].extension- rightFragFlapLength );
				  saveLocalAlignerVariables();
				  gotNewRightMA = GetNewUnitigMultiAlign( unitig, fragPoss, rFragIid );
				  restoreLocalAlignerVariables();
				  free( fragPoss );  // inefficient, let's just do a big array once that get's reused

				  if ( !gotNewRightMA )
					unitigMultiAlignFailures++;

				  if ( gotNewRightMA )
				  {
					int extendToLeft;
					SynchUnitigTWithMultiAlignT( unitig );
					
					fprintf( stderr, "before ReplaceEndUnitigInContig:\n");
					DumpContigMultiAlignInfo ( rcontig->id );
					  
					if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
					  extendToLeft = TRUE;
					else
					  extendToLeft = FALSE;

					saveLocalAlignerVariables();
					new_cma = ReplaceEndUnitigInContig( ScaffoldGraph->sequenceDB,
														ScaffoldGraph->fragStore,
														rcontig->id, unitig->id, extendToLeft, GlobalData->aligner);
					restoreLocalAlignerVariables();

					if ( !new_cma )
					  replaceEndUnitigFailures++;

					fprintf( stderr, "after ReplaceEndUnitigInContig:\n");
					DumpContigMultiAlignInfo ( rcontig->id );
					DumpContigMultiAlignInfoDirect ( new_cma, rcontig->id );
					  
					if ( new_cma )
					{
					  UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, rcontig->id, FALSE);
					  InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, rcontig->id, 
													 FALSE, new_cma, TRUE);
					}
					else
					  gotNewRightMA = FALSE;
					  
					fprintf( stderr, "after updating store:\n");
					DumpContigMultiAlignInfo ( rcontig->id );
					  
					fprintf( stderr, "strlen( Getchar (new_cma->consensus)): " F_SIZE_T "\n",
							 strlen( Getchar (new_cma->consensus, 0)));
					  
					// updateIntUnitigPoss( rcontig );
					if ( rcontigOrientation == A_B )
					  extendContig( rcontig, TRUE);
					else
					  extendContig( rcontig, FALSE);
				  }
				}
				  
				// to do: unwind changes if one or the other fails
				if (gotNewLeftMA == FALSE || gotNewRightMA == FALSE)  
				  keepGap = FALSE;
			  }
				
			  if (keepGap)
			  {
				double delta;
				int success;

				// now shift the right contig into place
				rcontig->offsetAEnd.mean += closedGapDelta[gapNumber];
				rcontig->offsetBEnd.mean += closedGapDelta[gapNumber];				  

				if (1)
				{
				  if (lcontigOrientation == A_B)
					fprintf( stderr, " after altering, lctg: %12.0f, %12.0f\n",
							 lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean);
				  else
					fprintf( stderr, " after altering, lctg: %12.0f, %12.0f\n",
							 lcontig->offsetBEnd.mean, lcontig->offsetAEnd.mean);
				  if (rcontigOrientation == A_B)
					fprintf( stderr, " after altering, rctg: %12.0f, %12.0f\n",
							 rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean);
				  else
					fprintf( stderr, " after altering, rctg: %12.0f, %12.0f\n",
							 rcontig->offsetBEnd.mean, rcontig->offsetAEnd.mean);
				}

/* temp hack!!*/  // if (lFragIid == -1 && rFragIid == -1)
				// this effectively kills the contigs and rebuilds from the unitigs
					
				// setup for contig merge
				if ( ContigPositions == NULL )
				  ContigPositions = CreateVA_IntElementPos(2);
				ResetVA_IntElementPos( ContigPositions );
				  
				if ( lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean )
				  delta = lcontig->offsetAEnd.mean;
				else
				  delta = lcontig->offsetBEnd.mean;
				  
				contigPos.ident = lcontig->id;
				contigPos.type = AS_CONTIG;
				contigPos.position.bgn = lcontig->offsetAEnd.mean - delta;
				contigPos.position.end = lcontig->offsetBEnd.mean - delta;
				AppendIntElementPos(ContigPositions, &contigPos);
				  
				fprintf( stderr, "lcontig %8d positioned at %8d, %8d\n", 
						 lcontig->id,contigPos.position.bgn, contigPos.position.end);
				  
				contigPos.ident = rcontig->id;
				contigPos.type = AS_CONTIG;
				contigPos.position.bgn = rcontig->offsetAEnd.mean - delta;
				contigPos.position.end = rcontig->offsetBEnd.mean - delta;
				AppendIntElementPos(ContigPositions, &contigPos);
				  
				fprintf( stderr, "rcontig %8d positioned at %8d, %8d\n", 
						 rcontig->id,contigPos.position.bgn, contigPos.position.end);

				if ( lcontigOrientation == A_B )
				  newOffsetAEnd.mean = lcontig->offsetAEnd.mean;
				else
				  newOffsetAEnd.mean = lcontig->offsetBEnd.mean;
				  
				if ( rcontigOrientation == A_B )
				  newOffsetBEnd.mean = rcontig->offsetBEnd.mean;
				else
				  newOffsetBEnd.mean = rcontig->offsetAEnd.mean;

#if 0
				{
				  MaxGaps = 5;
				  MaxBegGap = 200;
				  MaxEndGap = 200;
				  MaxInteriorGap = 30;
				  asymmetricEnds = TRUE;
				}
#endif				  
				fprintf( stderr, "before CreateAContigInScaffold:\n");
				DumpContigMultiAlignInfo ( lcontig->id );
				DumpContigMultiAlignInfo ( rcontig->id );
				// DumpContigUngappedOffsets( rcontig->id );
				  
				  
				// have to call this routine with normalized positions
				success = CreateAContigInScaffold( scaff, ContigPositions, newOffsetAEnd, newOffsetBEnd);
				  
				if ( !success )
				{
				  // why does this happen?
				  // now we have to unwind everything we've done up above???

				  keepGap = FALSE;
				  fprintf( stderr, "overlap found in extendClearRanges but not by CNS!\n");

				  createAContigFailures++;
				}
				  
	            rcontig = GetGraphNode( ScaffoldGraph->ContigGraph, rcontigID);
	            lcontig = GetGraphNode( ScaffoldGraph->ContigGraph, lcontigID);
				if ( keepGap )
				{
				  fprintf( stderr, "closed gap %8d, contigs %8d and %8d, fragIids %9d and %9d\n",
						   gapNumber, lcontig->id, rcontig->id, lFragIid, rFragIid);
					
				  closedGap[ gapNumber ] = TRUE;
				  // closedGapAhang[ gapNumber ] = ahang;
				  // closedGapOlapLength[ gapNumber ] = currLength;
				  // closedGapBhang[ gapNumber ] = bhang;
				  // closedGapLcontigBasesIntact[ gapNumber ] = lcontigBasesIntact;
				  // closedGapRcontigBasesIntact[ gapNumber ] = rcontigBasesIntact;
					
				  numGapsClosed++;
				  totalOlapVariance += currLength * currLength;
				  totalOlapLength += currLength;
				  totalOlapDiffs += currDiffs;
					
				  if (gapSize.mean < 100.0)
				  {
					numSmallGapsClosed++;
					numSmallGapsClosedThisScaff++;
				  }
				  else
				  {
					numLargeGapsClosed++;
					numLargeGapsClosedThisScaff++;
					if (gapSize.mean > maxGapSizeClosed)
					{
					  maxGapSizeClosed = gapSize.mean;
					  maxGapSizeClosedNumber = gapNumber;
					}
				  }
					
				  // newContig takes the place of what was the right contig
				  rcontig = newContig = GetGraphNode( ScaffoldGraph->ContigGraph, 
													  GetNumGraphNodes( ScaffoldGraph->ContigGraph ) - 1);
				  // adjustUnitigCoords( newContig );
					
				  // now we need to adjust contigs past the new contig, if not on end
				  if ( nextContigIndex != NULLINDEX )
				  {
					LengthT scaffoldDelta;
					  
					scaffoldDelta.mean = max( newContig->offsetAEnd.mean, newContig->offsetBEnd.mean) 
					  - maxRContigOffset;
					scaffoldDelta.variance = ComputeFudgeVariance( scaffoldDelta.mean );
					AddDeltaToScaffoldOffsets( ScaffoldGraph, scaff->id, nextContigIndex,
											   TRUE, FALSE, scaffoldDelta);
				  }
				}
			  }

			  // if we didn't close gap for whatever reason undo all the frag and unitig changes
			  if ( keepGap == FALSE)
			  {
				fprintf( stderr, "did not close gap %8d, contigs %8d and %8d\n",
						 gapNumber, lcontig->id, rcontig->id);
				closedGap[ gapNumber ] = FALSE;
				revertToCnsClearRange( lFragIid );
				revertToCnsClearRange( rFragIid );				  
				restoreFragAndUnitigData( lFragIid, rFragIid );
			  }
			}
			else
			{
			  fprintf( stderr, "did not close gap %8d, contigs %8d and %8d\n",
					   gapNumber, lcontig->id, rcontig->id);
			  closedGap[ gapNumber ] = FALSE;
			  noOverlapFound++;
			}
		  }
		}
	  }
	  fprintf( stderr, "after gapNumber %d:\n", gapNumber);
	  if (numSmallGaps > 0)
		fprintf( stderr, "             numSmallGaps: %d (closed %d, %.2f%%)\n", 
				 numSmallGaps, numSmallGapsClosed, 100.0 * numSmallGapsClosed / numSmallGaps);
	  if (numLargeGaps > 0)
		fprintf( stderr, "             numLargeGaps: %d (closed %d, %.2f%%)\n", 
				 numLargeGaps, numLargeGapsClosed, 100.0 * numLargeGapsClosed / numLargeGaps);
	  fprintf( stderr, "                  allGaps: %d (closed %d, %.2f%%)\n", 
			   numSmallGaps + numLargeGaps, numSmallGapsClosed + numLargeGapsClosed, 
			   100.0 * (numSmallGapsClosed + numLargeGapsClosed) / (numSmallGaps + numLargeGaps));
	  gapNumber++;
	  lcontig = rcontig;
	  lcontigID = lcontig->id;
	  rcontigID = lcontig->BEndNext;

	  fprintf( stderr, "at bottom of loop: lcontig->BEndNext: %d\n", lcontig->BEndNext);
	}
	
	// checkpointing
	sumScaffoldLengths += scaff->bpLength.mean;

	fprintf( stderr, "sumScaffoldLengths: %f, sumScaffoldLengthsLastCkp: %f\n",
			 sumScaffoldLengths, sumScaffoldLengthsLastCkp);
	fprintf( stderr, "scaffold stats, scaff %10d, smallGaps %8d closed %8d, largeGaps %8d closed %8d\n",
			 scaff->id, numSmallGapsThisScaff, numSmallGapsClosedThisScaff,
			 numLargeGapsThisScaff, numLargeGapsClosedThisScaff);

	numSmallGapsInScaffold[sid] = numSmallGapsThisScaff;
	numSmallGapsClosedInScaffold[sid] = numSmallGapsClosedThisScaff;
	numLargeGapsInScaffold[sid] = numLargeGapsThisScaff;
	numLargeGapsClosedInScaffold[sid] = numLargeGapsClosedThisScaff;
	numGapsInScaffold[sid] = numSmallGapsInScaffold[sid] + numLargeGapsInScaffold[sid];
	numGapsClosedInScaffold[sid] = numSmallGapsClosedInScaffold[sid] + numLargeGapsClosedInScaffold[sid];	


	if(numSmallGapsClosedInScaffold[sid]+numLargeGapsClosedInScaffold[sid]>0){
	  int status = RECOMPUTE_SINGULAR;
	  int recomputeIteration = 0;
	  while(recomputeIteration < 3 &&
		(status == RECOMPUTE_SINGULAR ||
		 status == RECOMPUTE_CONTIGGED_CONTAINMENTS))
	    {
	      // need to make sure scaffold is connected with trusted raw edges
	      MarkInternalEdgeStatus(ScaffoldGraph,
				     GetGraphNode(ScaffoldGraph->ScaffoldGraph,
						  sid),
				     PAIRWISECHI2THRESHOLD_CGW,
				     1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
				     TRUE, TRUE, 0, TRUE);

	      assert(IsScaffoldInternallyConnected(ScaffoldGraph,
						   GetGraphNode(ScaffoldGraph->ScaffoldGraph,
								sid),
						   ALL_EDGES));
	      
	      status =
		RecomputeOffsetsInScaffold(ScaffoldGraph,
					   GetGraphNode(ScaffoldGraph->ScaffoldGraph,
							sid),
					   TRUE, TRUE, FALSE);
	    }
	}

	if ( sumScaffoldLengths - sumScaffoldLengthsLastCkp > 90000000 )
	{
	  sumScaffoldLengthsLastCkp = sumScaffoldLengths;
	  fprintf( stderr, "checkpoint %d written during extendClearRanges, sumScaffoldLengths: %f\n",
			   ScaffoldGraph->checkPointIteration, sumScaffoldLengths);
	  fprintf( GlobalData->timefp, "checkpoint %d written during extendClearRanges, sumScaffoldLengths: %f\n",
			   ScaffoldGraph->checkPointIteration, sumScaffoldLengths);
	  CheckpointScaffoldGraph(ScaffoldGraph);
	  {
		char temp_buf[1024];
		int sysReturn;
		sprintf( temp_buf, "cp %s/db.frg %s/db.frg.%d", 
				 GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
		sysReturn = system( temp_buf );
		if ( sysReturn != -1)
		  fprintf( stderr, "copied frgStore to %s/db.frg.%d\n",
				   GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
		else
		  fprintf( stderr, "error encountered copying frgStore to %s/db.frg.%d\n",
				   GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration);		  

		writeEcrCheckpoint( numGapsInScaffold, numGapsClosedInScaffold,
							numSmallGapsInScaffold, numSmallGapsClosedInScaffold,
							numLargeGapsInScaffold, numLargeGapsClosedInScaffold);
	  }
	}
  }

  // Variance = mean(x^2) - (mean(x))^2
  if ( numGapsClosed > 0 )
	totalOlapVariance = (totalOlapVariance / numGapsClosed ) - 
	  ((float) totalOlapLength / numGapsClosed) * ((float) totalOlapLength / numGapsClosed);
  else
	totalOlapVariance = -0.0;

  fprintf( stderr, "\n");
  fprintf( stderr, "                  numGaps: %d\n", numGaps);
  fprintf( stderr, "        numExtendableGaps: %d (left: %d, right: %d, both: %d)\n", 
		   numExtendableGaps, leftContigExtendable, rightContigExtendable, bothContigsExtendable );
  fprintf( stderr, "            numGapsClosed: %d\n", numGapsClosed);
  fprintf( stderr, "             numSmallGaps: %d (closed %d, %.2f%%)\n", 
		   numSmallGaps, numSmallGapsClosed, 100.0 * numSmallGapsClosed / numSmallGaps);
  fprintf( stderr, "             numLargeGaps: %d (closed %d, %.2f%%)\n", 
		   numLargeGaps, numLargeGapsClosed, 100.0 * numLargeGapsClosed / numLargeGaps);
  fprintf( stderr, "                  allGaps: %d (closed %d, %.2f%%)\n", 
		   numSmallGaps + numLargeGaps, numSmallGapsClosed + numLargeGapsClosed, 
		   100.0 * (numSmallGapsClosed + numLargeGapsClosed) / (numSmallGaps + numLargeGaps));
  fprintf( stderr, "         maxGapSizeClosed: %.2f (gap number %d) \n", maxGapSizeClosed, maxGapSizeClosedNumber);
  fprintf( stderr, "       numGapsVarTooSmall: %d\n", numGapsVarTooSmall);
  fprintf( stderr, "        numClosingsTried : %d\n", numClosingsTried);
  fprintf( stderr, "unitigMultiAlignFailures : %d\n", unitigMultiAlignFailures);
  fprintf( stderr, "replaceEndUnitigFailures : %d\n", replaceEndUnitigFailures);
  fprintf( stderr, "    createAContigFailures: %d\n", createAContigFailures);
  fprintf( stderr, "           noOverlapFound: %d\n", noOverlapFound);

  if ( numGapsClosed > 0 )
	fprintf( stderr, "            avgOlapLength: %.2f\n", (float) totalOlapLength / numGapsClosed);
  fprintf( stderr, "         stddevOlapLength: %.2f\n", (float) sqrt( totalOlapVariance ));
  if ( numGapsClosed > 0 )
	fprintf( stderr, "             avgOlapDiffs: %.2f\n", (float) totalOlapDiffs / numGapsClosed);
  fprintf( stderr, "          surrogatesOnEnd: %d\n", surrogatesOnEnd);
  fprintf( stderr, "\n");
  fprintf( stderr, "                  MaxGaps: %d\n", MaxGaps);
  fprintf( stderr, "           MaxInteriorGap: %d\n", MaxInteriorGap);
  fprintf( stderr, "             CONTIG_BASES: %d\n", CONTIG_BASES);
  fprintf( stderr, "   totalContigsBaseChange: %d\n", totalContigsBaseChange);
  fprintf( stderr, "   totalBasesInClosedGaps: %d\n", totalBasesInClosedGaps);
  // scaffold base change should be the diff between originalsScaffolds.cgm and alteredScaffolds.cgm
  fprintf( stderr, "     scaffold base change: %d\n", totalContigsBaseChange - totalBasesInClosedGaps);
  fprintf( stderr, "\n");

  fprintf( GlobalData->timefp, "checkpoint %d written at end of extendClearRanges, sumScaffoldLengths: %f",
		   ScaffoldGraph->checkPointIteration, sumScaffoldLengths);
  CheckpointScaffoldGraph(ScaffoldGraph);
  {
	char temp_buf[1024];
	int sysReturn;
	sprintf( temp_buf, "cp %s/db.frg %s/db.frg.ecr.%d", 
			 GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
	sysReturn = system( temp_buf );
	if ( sysReturn != -1)
	  fprintf( stderr, "copied frgStore to %s/db.frg.%d\n",
			   GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
	else
	  fprintf( stderr, "error encountered copying frgStore to %s/db.frg.%d\n",
			   GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration);		  

	writeEcrCheckpoint( numGapsInScaffold, numGapsClosedInScaffold,
						numSmallGapsInScaffold, numSmallGapsClosedInScaffold,
						numLargeGapsInScaffold, numLargeGapsClosedInScaffold);
  }
  
#if 0
  // numDeletedContigs = alterContigs_old( numGaps, closedGap, closedGapDelta, 
  //						   lcontigIdGap, rcontigIdGap, lcontigLength, rcontigLength,
  //							   contigValid, allContigLengths, alteredScaffoldLengths);

  numDeletedContigs = alterContigs( numGaps, closedGap, closedGapDelta, 
								   lcontigIdGap, rcontigIdGap, lcontigLength, rcontigLength,
									closedGapAhang,
									closedGapOlapLength,
									closedGapBhang,
									closedGapLcontigBasesIntact,
									closedGapRcontigBasesIntact,
									contigValid, allContigLengths, alteredScaffoldLengths);

  // produceContigStatsConventional just produces the same stats as cgw
  produceContigStatsConventional();

  //produceContigStats( numGaps, closedGap, closedGapDelta, lcontigIdGap, rcontigIdGap, lcontigLength, rcontigLength,
  //				  contigValid, allContigLengths, numDeletedContigs);

  produceGapStats( numGaps, closedGap, closedGapDelta, originalGaps);

  produceScaffoldStats( alteredScaffoldLengths );
#endif

  t1 = time(0);
  fprintf( stderr, "====> Ending at %s\n", ctime(&t1));

  exit(0);
}

// since we mucked with the unitigs multialignment, reset the offsets and bpLengths of all the unitigs in the contig
// the ones in the ScaffoldGraph are not valid anymore
void adjustUnitigCoords( NodeCGW_T *contig )
{
  MultiAlignT *ma;
  int i;
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
  
  for ( i = 0; i < GetNumIntUnitigPoss( ma->u_list ); i++)
  {
	IntUnitigPos *pos = GetIntUnitigPos( ma->u_list, i);
	NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, pos->ident);
	MultiAlignT *uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE);

	fprintf( stderr, "before unitig %8d, bgn: %10d, end: %10d, length: %10d\n", 
			 unitig->id, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
	
	fprintf( stderr, "in adjustUnitigCoords, for unitig %d strlen( ma->consensus ) = " F_SIZE_T "\n",
			 unitig->id, strlen( Getchar( uma->consensus, 0) ));

	unitig->bpLength.mean = pos->position.end - pos->position.bgn;      // set length
	if ( unitig->offsetAEnd.mean < unitig->offsetBEnd.mean )          // ordering info is okay
	{
	  unitig->offsetAEnd.mean = pos->position.bgn;
	  unitig->offsetBEnd.mean = pos->position.end;
	}
	else 
	{
	  unitig->offsetAEnd.mean = pos->position.end;
	  unitig->offsetBEnd.mean = pos->position.bgn;
	}

	fprintf( stderr, " after unitig %8d, bgn: %10d, end: %10d, length: %10d\n", 
			 unitig->id, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
  }
}

int alterContigs( int numGaps, int *closedGap, int *closedGapDelta,
				  int *lcontigIdGap, int *rcontigIdGap, int *lcontigLength, int *rcontigLength,
				  int *closedGapAhang,
				  int *closedGapOlapLength,
				  int *closedGapBhang,
				  int *closedGapLcontigBasesIntact,
				  int *closedGapRcontigBasesIntact,
				  int *contigValid, int *allContigLengths, int *alteredScaffoldLengths)
{
  int i;
  int currentContigLength = 0, rcontigLastClosedGap = -1;
  int startContig = NULLINDEX, endContig = NULLINDEX;
  int previousGapClosed = FALSE;
  int numDeletedContigs = 0;
  LengthT gapSize;
  int numContigsMerged = 0;
  int sumMergedContigLengths = 0;
  
  for ( i = 0; i < numGaps; i++)
  {
	if ( closedGap[i] == TRUE )
	{
	  fprintf( stderr, "gapNumber: %d, lcontig: %6d (len: %6d), rcontig: %6d (len: %6d), closedGapDelta: %5d\n", 
			   i, lcontigIdGap[i], lcontigLength[i], rcontigIdGap[i], rcontigLength[i], closedGapDelta[i]);
	}
  }
  
  for ( i = 0; i < numGaps; i++)
  {
	NodeCGW_T *lcontig, *rcontig;
	
	if ( closedGap[i] == TRUE )
	{
	  int scaffoldID;

	  lcontig = GetGraphNode( ScaffoldGraph->ContigGraph, lcontigIdGap[i]);
	  rcontig = GetGraphNode( ScaffoldGraph->ContigGraph, rcontigIdGap[i]);
	  scaffoldID = lcontig->scaffoldID;
	  
	  gapSize = FindGapLength( lcontig, rcontig, FALSE);

	  contigValid[ lcontigIdGap[i] ] = FALSE;   // mark that the left contig is no longer valid
	  numDeletedContigs++;
	  // alteredScaffoldLengths[ scaffoldID ] += (int) gapSize.mean + closedGapDelta[ i ];
	  
	  if ( lcontigIdGap[i] == rcontigLastClosedGap )
	  {
		// the "- CONTIG_BASES" term accounts for the fact that some of the bases of what is currently
		// the lcontig were counted already when it was part of the right contig.
		currentContigLength += - CONTIG_BASES + closedGapAhang[i] + closedGapOlapLength[i] +
		  closedGapBhang[i] + closedGapRcontigBasesIntact[i];
		rcontigLastClosedGap = rcontigIdGap[ i ];
		endContig = rcontigIdGap[ i ];
		numContigsMerged++;
		sumMergedContigLengths += (int) rcontig->bpLength.mean;
	  }
	  else
	  {
		if ( previousGapClosed )
		{
		  fprintf( stderr, "currentContigLength: %d (start contig: %d, end contig: %d, %d contigs merged)\n",
				   currentContigLength, startContig, endContig, numContigsMerged);
		  fprintf( stderr, "currentContigLength: %d, sumMergedContigLengths: %d, diff: %d\n", 
				   currentContigLength, sumMergedContigLengths, currentContigLength - sumMergedContigLengths);
		  allContigLengths[ endContig ] = currentContigLength;  // set the endContig to hold the length of merged contigs
		}
		// start anew
		currentContigLength = closedGapLcontigBasesIntact[i] + closedGapAhang[i] + closedGapOlapLength[i] +
		  closedGapBhang[i] + closedGapRcontigBasesIntact[i];
		rcontigLastClosedGap = rcontigIdGap[ i ];
		startContig = lcontigIdGap[ i ];
		endContig = rcontigIdGap[ i ];
		numContigsMerged = 2;
		sumMergedContigLengths = (int) lcontig->bpLength.mean + (int) rcontig->bpLength.mean;
	  }
	  previousGapClosed = TRUE;
	}
	else
	{
	  if ( previousGapClosed == TRUE )
	  {
		fprintf( stderr, "currentContigLength: %d (start contig: %d, end contig: %d, %d contigs merged)\n",
				 currentContigLength, startContig, endContig, numContigsMerged);
		fprintf( stderr, "currentContigLength: %d, sumMergedContigLengths: %d, diff: %d\n", 
				 currentContigLength, sumMergedContigLengths, currentContigLength - sumMergedContigLengths);
		allContigLengths[ endContig ] = currentContigLength;  // set the endContig to hold the length of merged contigs
	  }
	  rcontigLastClosedGap = -1;	  	
	  previousGapClosed = FALSE;
	  numContigsMerged = 0;
	}
  }
  return numDeletedContigs;
}

// int alterContigs_old( int numGaps, int *closedGap, int *closedGapDelta,
// 				  int *lcontigIdGap, int *rcontigIdGap, int *lcontigLength, int *rcontigLength,
// 				  int *contigValid, int *allContigLengths, int *alteredScaffoldLengths)
// {
//   int i;
//   int currentContigLength, rcontigLastClosedGap = -1;
//   int startContig, endContig;
//   int previousGapClosed = FALSE;
//   int numDeletedContigs = 0;
//   LengthT gapSize;
//   
//   for ( i = 0; i < numGaps; i++)
//   {
// 	if ( closedGap[i] == TRUE )
// 	{
// 	  fprintf( stderr, "gapNumber: %d, lcontig: %6d (len: %6d), rcontig: %6d (len: %6d), closedGapDelta: %5d\n", 
// 			   i, lcontigIdGap[i], lcontigLength[i], rcontigIdGap[i], rcontigLength[i], closedGapDelta[i]);
// 	}
//   }
//   
//   for ( i = 0; i < numGaps; i++)
//   {
// 	NodeCGW_T *lcontig, *rcontig;
// 	
// 	if ( closedGap[i] == TRUE )
// 	{
// 	  int scaffoldID;
// 
// 	  lcontig = GetGraphNode( ScaffoldGraph->ContigGraph, lcontigIdGap[i]);
// 	  rcontig = GetGraphNode( ScaffoldGraph->ContigGraph, rcontigIdGap[i]);
// 	  scaffoldID = lcontig->scaffoldID;
// 	  
// 	  gapSize = FindGapLength( lcontig, rcontig, FALSE);
// 
// 	  contigValid[ lcontigIdGap[i] ] = FALSE;   // mark that the left contig is no longer valid
// 	  numDeletedContigs++;
// 	  alteredScaffoldLengths[ scaffoldID ] += (int) gapSize.mean + closedGapDelta[ i ];
// 	  
// 	  if ( lcontigIdGap[i] == rcontigLastClosedGap )
// 	  {
// 		currentContigLength += (int) gapSize.mean + closedGapDelta[ i ] + rcontigLength[ i ];
// 		endContig = rcontigIdGap[ i ];
// 	  }
// 	  else
// 	  {
// 		if ( previousGapClosed )
// 		{
// 		  fprintf( stderr, "currentContigLength: %d (start contig: %d, end contig: %d)\n",
// 				   currentContigLength, startContig, endContig);
// 		  allContigLengths[ endContig ] = currentContigLength;  // set the endContig to hold the length of merged contigs
// 		}
// 		// start anew
// 		currentContigLength = lcontigLength[ i ] + (int) gapSize.mean + closedGapDelta[ i ] + rcontigLength[ i ];
// 		rcontigLastClosedGap = rcontigIdGap[ i ];
// 		startContig = lcontigIdGap[ i ];
// 		endContig = rcontigIdGap[ i ];
// 	  }
// 	  previousGapClosed = TRUE;
// 	}
// 	else
// 	{
// 	  if ( previousGapClosed == TRUE )
// 	  {
// 		fprintf( stderr, "currentContigLength: %d (start contig: %d, end contig: %d)\n",
// 				 currentContigLength, startContig, endContig);
// 		allContigLengths[ endContig ] = currentContigLength;  // set the endContig to hold the length of merged contigs
// 	  }
// 	  rcontigLastClosedGap = -1;	  	
// 	  previousGapClosed = FALSE;
// 	}
//   }
//   return numDeletedContigs;
// }

void collectContigStats( int *allContigLengths, int *contigValid)
{
  double totalLengthContigs = 0.0;
  int numContigs = 0;
  int sid;
  
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
	CIScaffoldTIterator CIsTemp;
	CIScaffoldT * scaff;
	int icnt;
	// int leftFragFound, rightFragFound;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    
	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    
	if ((isDeadCIScaffoldT(scaff)) ||
		(scaff->type != REAL_SCAFFOLD)) // || (scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}
	fprintf(stderr,"\n=====================================================================\n");
	fprintf(stderr,"=== examing scaffold %d, size %f\n", sid, scaff->bpLength.mean);

	// make sure the scaffold is there
	assert(scaff != NULL);

	icnt = 0;
	InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,	FALSE, &CIsTemp);
	while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ContigT *contig;
	  
	  contig = GetGraphNode( ScaffoldGraph->RezGraph, CIsTemp.curr);
	  assert(contig != NULL);

	  contigValid[ contig->id ] = TRUE;
	  allContigLengths[ contig->id ] = (int) contig->bpLength.mean;

	  totalLengthContigs += contig->bpLength.mean;
	  numContigs++;

	  // if at end of scaffold break
	  if (CIsTemp.next == -1)
		break;
	}
  }
  fprintf( stderr, "collectContigStats, numContigs: %d, avg contig length: %f\n", 
		   numContigs, totalLengthContigs / numContigs);
}

void produceContigStatsConventional()
{
  double totalLengthContigs = 0.0;
  int numContigs = 0;
  int sid;
  
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
	CIScaffoldTIterator CIsTemp;
	CIScaffoldT * scaff;
	int icnt;
	// int leftFragFound, rightFragFound;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    
	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    
	if ((isDeadCIScaffoldT(scaff)) ||
		(scaff->type != REAL_SCAFFOLD)) // || (scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}
	fprintf(stderr,"\n=====================================================================\n");
	fprintf(stderr,"=== examing scaffold %d, size %f\n", sid, scaff->bpLength.mean);

	// make sure the scaffold is there
	assert(scaff != NULL);

	icnt = 0;
	InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,	FALSE, &CIsTemp);
	while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ContigT *contig;
	  
	  contig = GetGraphNode( ScaffoldGraph->RezGraph, CIsTemp.curr);
	  assert(contig != NULL);

	  totalLengthContigs += contig->bpLength.mean;
	  numContigs++;

	  // if at end of scaffold break
	  if (CIsTemp.next == -1)
		break;
	}
  }
  fprintf( stderr, "produceContigStatsConventional, numContigs: %d, avg contig length: %f\n", 
		   numContigs, totalLengthContigs / numContigs);
}

void produceGapStats( int numGaps, int *closedGap, int *closedGapDelta, int *originalGaps)
{
  int i;
  FILE *alteredGapsFile;
  FILE *originalGapsFile;
  
  originalGapsFile = fopen( "originalGaps.cgm", "w");
  if ( originalGapsFile == NULL )
  {
	fprintf( stderr, "failed to open originalGapsFile\n");
	assert(0);
  }
  fprintf( originalGapsFile, "original gaps\n");

  for ( i = 0; i < numGaps; i++)
	fprintf( originalGapsFile, "%d\n", originalGaps[i]);

  fclose( originalGapsFile );

  alteredGapsFile = fopen( "alteredGaps.cgm", "w");
  if ( alteredGapsFile == NULL )
  {
	fprintf( stderr, "failed to open alteredGapsFile\n");
	assert(0);
  }
  fprintf( alteredGapsFile, "altered gaps\n");

  for ( i = 0; i < numGaps; i++)
  {
	if ( closedGap[i] == FALSE)
	  fprintf( alteredGapsFile, "%d\n", originalGaps[i]);
  }
  
  fclose( originalGapsFile );
}

void produceScaffoldStats( int *alteredScaffoldLengths )
{
  int sid;
  FILE *alteredScaffoldsFile;
  FILE *originalScaffoldsFile;
  
  originalScaffoldsFile = fopen( "originalScaffolds.cgm", "w");
  if ( originalScaffoldsFile == NULL )
  {
	fprintf( stderr, "failed to open originalScaffoldsFile\n");
	assert(0);
  }
  fprintf( originalScaffoldsFile, "original scaffolds\n");

  alteredScaffoldsFile = fopen( "alteredScaffolds.cgm", "w");
  if ( alteredScaffoldsFile == NULL )
  {
	fprintf( stderr, "failed to open alteredScaffoldsFile\n");
	assert(0);
  }
  fprintf( alteredScaffoldsFile, "altered scaffolds\n");

  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
	CIScaffoldT * scaff;
	// int leftFragFound, rightFragFound;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    
	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    
	if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD) 
		|| (scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}
	
	// make sure the scaffold is there
	assert(scaff != NULL);
	
	fprintf( originalScaffoldsFile, "%d\n", (int) scaff->bpLength.mean);
	// at this point alteredScaffoldLengths has only the deltas for its gaps, not the deltas plus the original length 
	fprintf( alteredScaffoldsFile, "%d\n", alteredScaffoldLengths[sid]);
  }
  
  fclose( originalScaffoldsFile );
  fclose( alteredScaffoldsFile );
}

void produceContigStats( int numGaps, int *closedGap, int *closedGapDelta,
                         int *lcontigIdGap, int *rcontigIdGap, int *lcontigLength, int *rcontigLength,
                         int *contigValid, int *allContigLengths, int numDeletedContigs)
{
  double totalLengthContigs = 0.0;
  int numContigs = 0;
  int sid;
  FILE *alteredContigsFile;
  
  alteredContigsFile = fopen( "alteredContigs.cgm", "w");
  if ( alteredContigsFile == NULL )
  {
	fprintf( stderr, "failed to open alteredContigsFile\n");
	assert(0);
  }
  fprintf( alteredContigsFile, "altered contigs\n");

  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
	CIScaffoldTIterator CIsTemp;
	CIScaffoldT * scaff;
	int icnt;
	// int leftFragFound, rightFragFound;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    
	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    
	if ((isDeadCIScaffoldT(scaff)) ||
		(scaff->type != REAL_SCAFFOLD)) // || (scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}
	fprintf(stderr,"\n=====================================================================\n");
	fprintf(stderr,"=== examing scaffold %d, size %f\n", sid, scaff->bpLength.mean);

	// make sure the scaffold is there
	assert(scaff != NULL);

	icnt = 0;
	InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,	FALSE, &CIsTemp);
	while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ContigT *contig;
	  
	  contig = GetGraphNode( ScaffoldGraph->RezGraph, CIsTemp.curr);
	  assert(contig != NULL);

	  // totalLengthContigs += contig->bpLength.mean;
	  if ( contigValid[ contig->id ] )
	  {
		totalLengthContigs += allContigLengths[ contig->id ];
		numContigs++;
		fprintf( alteredContigsFile, "%d\n", allContigLengths[ contig->id ]);
	  }
	  
	  // if at end of scaffold break
	  if (CIsTemp.next == -1)
		break;
	}
  }
  fclose( alteredContigsFile );
  fprintf( stderr, "produceContigStats, numContigs: %d, avg contig length: %f\n", 
		   numContigs, totalLengthContigs / numContigs);
}

void initVarArrays(void)
{
  if ( lContigConsensus == NULL )
  {
	lContigConsensus = CreateVA_char(1024);
	rContigConsensus = CreateVA_char(1024);
	lContigQuality = CreateVA_char(1024);
	rContigQuality = CreateVA_char(1024);
	reformed_consensus = CreateVA_char(200000);
	reformed_quality = CreateVA_char(200000);
	reformed_deltas = CreateVA_int32(1);
  }
}

// findFirstFrag looks for a 3p->5p frag at the low end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int findFirstFrag( ContigT *contig, int *fragIid, int *extensionOut, int *basesToNextFrag )
{
  MultiAlignT *ma;
  IntMultiPos *mp;
  CIFragT *frag;
  int i, numFrags;
  int foundFrag, fragOnEnd = FALSE;
  int currExtension = 0, secondFragStart = 10000;
  
  fprintf( stderr, "in FindFirstFrag\n");
  
  ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss( ma->f_list);

  // fprintf( stderr, "contig %d has %d frags over a length of %f\n", contig->id, numFrags, contig->bpLength.mean);
  // fprintf( stderr, "contig %d has %d unitigs over a length of %f\n", contig->id, 
  //   GetNumIntUnitigPoss(ma->u_list), contig->bpLength.mean);
  
  foundFrag = FALSE;
  *extensionOut = 0;
  for ( i = 0; i < numFrags; i++)
  {
	mp = GetIntMultiPos( ma->f_list, i);
	frag = GetCIFragT( ScaffoldGraph->CIFrags, (int32) mp->source);

	if (frag->contigOffset3p.mean < 100.0 &&      // frag is within a cutoff of the low end of the contig
		frag->locale == -1 &&                     // and is a read
		frag->contigOffset3p.mean < frag->contigOffset5p.mean)  // and points in the right direction
	{
	  char seqbuffer[AS_BACTIG_MAX_LEN+1], qltbuffer[AS_BACTIG_MAX_LEN+1];
	  unsigned int clr_bgn, clr_end;
	  int frag3pExtra, extension;
	  
	  if ( fsread == NULL ) fsread = new_ReadStruct();

	  getFragStore( ScaffoldGraph->fragStore, frag->iid, FRAG_S_ALL, fsread);
	  getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	  getSequence_ReadStruct( fsread, seqbuffer, qltbuffer, AS_BACTIG_MAX_LEN);


	  //                 <--------------------------------------------------------------------------- contig
	  // 3p <------------------|---------------------------------------------|----- 5p frag
	  //                    clr_end                                       clr_bgn
	  //    |-----------|
	  //      extension

	  //             <--------------------------------------------------------------------------- contig
	  //                     <-|---------------------------------------------|----- 5p frag
	  //                    clr_end                                       clr_bgn
	  //             |------|
	  //               extension (negative)

	  fprintf( stderr, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
	  fprintf( stderr, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
			   frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
	  fprintf( stderr, "frag length: " F_SIZE_T ", 3p past clr_end length: " F_SIZE_T "\n", strlen( seqbuffer ), 
			   strlen( seqbuffer ) - clr_end);
	  fprintf( stderr, "extension: " F_SIZE_T "\n", strlen( seqbuffer ) - clr_end - (int) frag->contigOffset3p.mean);
	  
	  frag3pExtra = strlen( seqbuffer ) - clr_end;
	  extension = frag3pExtra - frag->contigOffset3p.mean;

	  if ( extension > currExtension )
	  {
		currExtension = strlen( seqbuffer ) - clr_end - frag->contigOffset3p.mean;
		foundFrag = TRUE;
		if ( frag->contigOffset3p.mean == 0 )
		  fragOnEnd = TRUE;
		else
		  fragOnEnd = FALSE;
		*extensionOut = currExtension;
		*fragIid = frag->iid;

		if (debug == 1)
		{
		  fprintf( stderr, "in contig %d, frag %d is at %f -> %f (5p->3p) \n", 
				   contig->id, frag->iid,
				   frag->contigOffset5p.mean, frag->contigOffset3p.mean);
		  fprintf( stderr, "extension ratio: %.2f\n", extension / (float) ( 1.0 + frag3pExtra - extension ));
		}
	  }
	}

	// secondFragStart is where the next to end frag starts, and thus where we start 2x coverage
	// we don't care if it's the 3p or 5p end
	if ( frag->contigOffset3p.mean > 0 && (int) frag->contigOffset3p.mean < secondFragStart)
	{
	  secondFragStart = frag->contigOffset3p.mean;
	  // fprintf( stderr, "secondFragStart %d set by frag %d (3p)\n", secondFragStart, frag->iid);
	}
	if ( (int) frag->contigOffset5p.mean < secondFragStart)
	{
	  secondFragStart = frag->contigOffset5p.mean;
	  // fprintf( stderr, "secondFragStart %d set by frag %d (5p)\n", secondFragStart, frag->iid);
	}
  }

  if ( fragOnEnd == TRUE )
	*basesToNextFrag = secondFragStart;
  else
	*basesToNextFrag = 0;

  // fprintf( stderr, "basesToNextFrag: %d\n", *basesToNextFrag);
  
  return foundFrag;
}

#if 1

// findFirstFrag looks for a 3p->5p frag at the low end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int findFirstExtendableFrags( ContigT *contig, extendableFrag *extFragsArray )
{
  MultiAlignT *ma;
  IntMultiPos *mp;
  CIFragT *frag;
  int i, numFrags, firstUnitigID;
  int extendableFragCount = 0;
  int secondFragStart = 10000;

  ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss( ma->f_list);
  firstUnitigID = GetIntUnitigPos( ma->u_list, 0 )->ident;

  fprintf( stderr, "in findFirstExtendableFrags, firstUnitigID: %d\n", firstUnitigID);
  
  // foundFrag = FALSE;
  // *extensionOut = 0;
  for ( i = 0; i < numFrags; i++)
  {
	mp = GetIntMultiPos( ma->f_list, i);
	frag = GetCIFragT( ScaffoldGraph->CIFrags, (int32) mp->source);

	if (frag->contigOffset3p.mean < 100.0 &&      // frag is within a cutoff of the low end of the contig
		frag->locale == -1 &&                     // and is a read
		frag->contigOffset3p.mean < frag->contigOffset5p.mean &&  // and points in the right direction
		frag->cid == firstUnitigID ) // and is in the first unitig
	{
	  char seqbuffer[AS_BACTIG_MAX_LEN+1], qltbuffer[AS_BACTIG_MAX_LEN+1];
	  unsigned int clr_bgn, clr_end;
	  int frag3pExtra, extension;
	  
	  if ( fsread == NULL ) fsread = new_ReadStruct();

	  getFragStore( ScaffoldGraph->fragStore, frag->iid, FRAG_S_ALL, fsread);
	  getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	  getSequence_ReadStruct( fsread, seqbuffer, qltbuffer, AS_BACTIG_MAX_LEN);


	  //                 <--------------------------------------------------------------------------- contig
	  // 3p <------------------|---------------------------------------------|----- 5p frag
	  //                    clr_end                                       clr_bgn
	  //    |-----------|
	  //      extension

	  //             <--------------------------------------------------------------------------- contig
	  //                     <-|---------------------------------------------|----- 5p frag
	  //                    clr_end                                       clr_bgn
	  //             |------|
	  //               extension (negative)

	  fprintf( stderr, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
	  fprintf( stderr, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
			   frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
	  fprintf( stderr, "frag length: " F_SIZE_T ", 3p past clr_end length: " F_SIZE_T "\n", strlen( seqbuffer ), 
			   strlen( seqbuffer ) - clr_end);
	  fprintf( stderr, "extension: " F_SIZE_T "\n", strlen( seqbuffer ) - clr_end - (int) frag->contigOffset3p.mean);
	  
	  frag3pExtra = strlen( seqbuffer ) - clr_end;
	  extension = frag3pExtra - frag->contigOffset3p.mean;

	  // ask Granger what min extension we should accept
	  if ( extension > 30 )
	  {
		// foundFrag = TRUE;

		extFragsArray[ extendableFragCount ].fragIid = frag->iid;
		extFragsArray[ extendableFragCount ].extension = extension;
		extFragsArray[ extendableFragCount ].addedBases = frag3pExtra;

		fprintf( stderr, "for frag %d, extension: %8d, frag3pExtra: %8d\n",
				 frag->iid, extension, frag3pExtra);

		if ( frag->contigOffset3p.mean == 0 )
		  extFragsArray[ extendableFragCount ].fragOnEnd = TRUE;
		else
		  extFragsArray[ extendableFragCount ].fragOnEnd = FALSE;

		if (debug == 1)
		{
		  fprintf( stderr, "in contig %d, frag %d is at %f -> %f (5p->3p) \n", 
				   contig->id, frag->iid,
				   frag->contigOffset5p.mean, frag->contigOffset3p.mean);
		  fprintf( stderr, "extension ratio: %.2f\n", extension / (float) ( 1.0 + frag3pExtra - extension ));
		}

		extendableFragCount++;
		if ( extendableFragCount > MAX_EXTENDABLE_FRAGS)
		{
		  fprintf( stderr, "extendableFragCount (%d) is greater than MAX_EXTENDABLE_FRAGS, aborting...\n",
				   extendableFragCount);
		  assert(0);
		}
	  }
	}

	// secondFragStart is where the next to end frag starts, and thus where we start 2x coverage
	// we don't care if it's the 3p or 5p end
	if ( frag->contigOffset3p.mean > 0 && (int) frag->contigOffset3p.mean < secondFragStart)
	{
	  secondFragStart = frag->contigOffset3p.mean;
	  // fprintf( stderr, "secondFragStart %d set by frag %d (3p)\n", secondFragStart, frag->iid);
	}
	if ( (int) frag->contigOffset5p.mean < secondFragStart)
	{
	  secondFragStart = frag->contigOffset5p.mean;
	  // fprintf( stderr, "secondFragStart %d set by frag %d (5p)\n", secondFragStart, frag->iid);
	}
  }
  
  // now sort the extendable frags by their extendability
  qsort( extFragsArray, extendableFragCount, sizeof( extendableFrag ), &compExtendableFrags );

  fprintf( stderr, "extendableFragCount: %d\n", extendableFragCount);
  for ( i = 0; i < extendableFragCount; i++)
  {
	if (extFragsArray[i].fragOnEnd == TRUE)
	  extFragsArray[i].basesToNextFrag = secondFragStart;
	else
	  extFragsArray[i].basesToNextFrag = 0;

	fprintf( stderr, "contig %8d, frag %8d can extend %8d bases into the gap\n",
			 contig->id, extFragsArray[i].fragIid, extFragsArray[i].extension);
  }

  return extendableFragCount;
}
#endif

// findLastFrag looks for a 5p->3p frag at the high end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int findLastFrag( ContigT *contig, int *fragIid, int *extensionOut, int *basesToNextFrag )
{
  MultiAlignT *ma;
  IntMultiPos *mp;
  CIFragT *frag;
  int i, numFrags;
  int foundFrag = FALSE, fragOnEnd = FALSE;
  float maxContigPos;
  int currExtension = 0, secondFragEnd = 0;

  fprintf( stderr, "in FindLastFrag\n");
  
  ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss( ma->f_list);
  
  maxContigPos = contig->bpLength.mean - 1.0;
  
  *extensionOut = 0;
  for ( i = 0; i < numFrags; i++)
  {
	mp = GetIntMultiPos( ma->f_list, i);
	frag = GetCIFragT( ScaffoldGraph->CIFrags, (int32) mp->source);

	// if (frag->contigOffset3p.mean == maxContigPos && frag->locale == -1)
	if (frag->contigOffset3p.mean > maxContigPos - 100.0 &&      // frag is within a cutoff of the high end of the contig
		frag->locale == -1 &&                                    // and is a read
		frag->contigOffset5p.mean < frag->contigOffset3p.mean)   // and points in the right direction
	{
	  char seqbuffer[AS_BACTIG_MAX_LEN+1], qltbuffer[AS_BACTIG_MAX_LEN+1];
	  unsigned int clr_bgn, clr_end;
	  int frag3pExtra, extension;

	  if ( fsread == NULL ) fsread = new_ReadStruct();

	  getFragStore( ScaffoldGraph->fragStore, frag->iid, FRAG_S_ALL, fsread);
	  getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	  getSequence_ReadStruct( fsread, seqbuffer, qltbuffer, AS_BACTIG_MAX_LEN);

	  //    contig ----------------------------------------------------------------------------------->
	  //                                  5p -------|---------------------------------------------|------------> 3p 
	  //                                         clr_bgn                                       clr_end
	  //                                                                                               |-------|
	  //                                                                                             extension
	  frag3pExtra = strlen( seqbuffer ) - clr_end;
	  extension = frag3pExtra - (int) (contig->bpLength.mean - frag->contigOffset3p.mean);

	  fprintf( stderr, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
	  fprintf( stderr, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
			   frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
	  fprintf( stderr, "frag length: " F_SIZE_T ", 3p past clr_end length: %d\n", strlen( seqbuffer ), frag3pExtra);
	  fprintf( stderr, "extension: %d\n", extension);
	  
	  if ( extension > currExtension )
	  {
		currExtension = extension;
		foundFrag = TRUE;
		if (frag->contigOffset3p.mean == contig->bpLength.mean)
		  fragOnEnd = TRUE;
		else
		  fragOnEnd = FALSE;
		*extensionOut = currExtension;
		*fragIid = frag->iid;
		if (debug == 1)
		{
		  fprintf( stderr, "in contig %d, frag %d is at %f -> %f (5p->3p) maxContigPos: %f\n", 
				   contig->id, frag->iid,
				   frag->contigOffset5p.mean, frag->contigOffset3p.mean, maxContigPos);
		  fprintf( stderr, "extension ratio: %.2f\n", extension / (float) ( 1.0 + frag3pExtra - extension ));
		}
	  }
	}

	// secondFragEnd is where the next to end frag ends, and thus where we have 2x coverage
	// we don't care if it's the 3p or 5p end
	if ( frag->contigOffset3p.mean < contig->bpLength.mean && (int) frag->contigOffset3p.mean > secondFragEnd )
	{
	  secondFragEnd = (int) frag->contigOffset3p.mean;
	  // fprintf( stderr, "secondFragEnd %d set by frag %d (3p)\n", secondFragEnd, frag->iid);
	}
	if ( (int) frag->contigOffset5p.mean > secondFragEnd)
	{
	  secondFragEnd = frag->contigOffset5p.mean;
	  // fprintf( stderr, "secondFragEnd %d set by frag %d (5p)\n", secondFragEnd, frag->iid);
	}
  }
  if ( fragOnEnd == TRUE )
	*basesToNextFrag = (int) contig->bpLength.mean - secondFragEnd;
  else
	*basesToNextFrag = 0;
  // fprintf( stderr, "basesToNextFrag: %d\n", *basesToNextFrag);

  return foundFrag;
}

// findLastFrag looks for a 5p->3p frag at the high end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int findLastExtendableFrags( ContigT *contig, extendableFrag *extFragsArray )
{
  MultiAlignT *ma;
  IntMultiPos *mp;
  CIFragT *frag;
  int i, numFrags, lastUnitigID;
  float maxContigPos;
  int secondFragEnd = 0;
  int extendableFragCount = 0;
  // int basesToNextFrag;
  
  ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss( ma->f_list);
  lastUnitigID = GetIntUnitigPos( ma->u_list, GetNumIntUnitigPoss(ma->u_list) - 1 )->ident;  
  maxContigPos = contig->bpLength.mean - 1.0;

  fprintf( stderr, "in FindLastExtendableFrags, lastUnitigID: %d\n", lastUnitigID);
  
  // *extensionOut = 0;
  for ( i = 0; i < numFrags; i++)
  {
	mp = GetIntMultiPos( ma->f_list, i);
	frag = GetCIFragT( ScaffoldGraph->CIFrags, (int32) mp->source);

	// if (frag->contigOffset3p.mean == maxContigPos && frag->locale == -1)
	if (frag->contigOffset3p.mean > maxContigPos - 100.0 &&      // frag is within a cutoff of the high end of the contig
		frag->locale == -1 &&                                    // and is a read
		frag->contigOffset5p.mean < frag->contigOffset3p.mean && // and points in the right direction
		frag->cid == lastUnitigID)                               // and is in the last unitig
	{
	  char seqbuffer[AS_BACTIG_MAX_LEN+1], qltbuffer[AS_BACTIG_MAX_LEN+1];
	  unsigned int clr_bgn, clr_end;
	  int frag3pExtra, extension;

	  if ( fsread == NULL ) fsread = new_ReadStruct();

	  getFragStore( ScaffoldGraph->fragStore, frag->iid, FRAG_S_ALL, fsread);
	  getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	  getSequence_ReadStruct( fsread, seqbuffer, qltbuffer, AS_BACTIG_MAX_LEN);

	  //    contig ----------------------------------------------------------------------------------->
	  //                                  5p -------|---------------------------------------------|------------> 3p 
	  //                                         clr_bgn                                       clr_end
	  //                                                                                               |-------|
	  //                                                                                             extension
	  frag3pExtra = strlen( seqbuffer ) - clr_end;
	  extension = frag3pExtra - (int) (contig->bpLength.mean - frag->contigOffset3p.mean);

	  fprintf( stderr, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
	  fprintf( stderr, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
			   frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
	  fprintf( stderr, "frag length: " F_SIZE_T ", 3p past clr_end length: %d\n", strlen( seqbuffer ), frag3pExtra);
	  fprintf( stderr, "extension: %d\n", extension);
	  
	  if ( extension > 30 )
	  {
		// foundFrag = TRUE;

		extFragsArray[ extendableFragCount ].fragIid = frag->iid;
		extFragsArray[ extendableFragCount ].extension = extension;
		extFragsArray[ extendableFragCount ].addedBases = frag3pExtra;

		fprintf( stderr, "for frag %d, extension: %8d, frag3pExtra: %8d\n",
				 frag->iid, extension, frag3pExtra);

		if (frag->contigOffset3p.mean == contig->bpLength.mean)
		  extFragsArray[ extendableFragCount ].fragOnEnd = TRUE;
		else
		  extFragsArray[ extendableFragCount ].fragOnEnd = FALSE;
		if (debug == 1)
		{
		  fprintf( stderr, "in contig %d, frag %d is at %f -> %f (5p->3p) maxContigPos: %f\n", 
				   contig->id, frag->iid,
				   frag->contigOffset5p.mean, frag->contigOffset3p.mean, maxContigPos);
		  fprintf( stderr, "extension ratio: %.2f\n", extension / (float) ( 1.0 + frag3pExtra - extension ));
		}

		extendableFragCount++;
		if ( extendableFragCount > MAX_EXTENDABLE_FRAGS)
		{
		  fprintf( stderr, "extendableFragCount (%d) is greater than MAX_EXTENDABLE_FRAGS, aborting...\n",
				   extendableFragCount);
		  assert(0);
		}
	  }
	}

	// secondFragEnd is where the next to end frag ends, and thus where we have 2x coverage
	// we don't care if it's the 3p or 5p end
	if ( frag->contigOffset3p.mean < contig->bpLength.mean && (int) frag->contigOffset3p.mean > secondFragEnd )
	{
	  secondFragEnd = (int) frag->contigOffset3p.mean;
	  // fprintf( stderr, "secondFragEnd %d set by frag %d (3p)\n", secondFragEnd, frag->iid);
	}
	if ( (int) frag->contigOffset5p.mean > secondFragEnd)
	{
	  secondFragEnd = frag->contigOffset5p.mean;
	  // fprintf( stderr, "secondFragEnd %d set by frag %d (5p)\n", secondFragEnd, frag->iid);
	}
  }

  // now sort the extendable frags by their extendability
  qsort( extFragsArray, extendableFragCount, sizeof( extendableFrag ), &compExtendableFrags );

  fprintf( stderr, "extendableFragCount: %d\n", extendableFragCount);
  for ( i = 0; i < extendableFragCount; i++)
  {
	if (extFragsArray[i].fragOnEnd == TRUE)
	  extFragsArray[i].basesToNextFrag = secondFragEnd;
	else
	  extFragsArray[i].basesToNextFrag = 0;	  
	
	fprintf( stderr, "contig %8d, frag %8d can extend %8d bases into the gap\n",
			 contig->id, extFragsArray[i].fragIid, extFragsArray[i].extension);
  }

  return extendableFragCount;
}

// findLastUnitig looks for a surrogate at the high end of a contig
int findLastUnitig( ContigT *contig, int *unitigID )
{
  MultiAlignT *ma;
  int i, numUnitigs, isSurrogate = FALSE;
  float maxContigPos = 0.0;
  NodeCGW_T *unitig = NULL;
  
  fprintf( stderr, "in FindLastUnitig\n");
  
  ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numUnitigs = GetNumIntUnitigPoss( ma->u_list );
  
  // maxContigPos = contig->bpLength.mean - 1.0;
  
  // can't just jump to last unitig since unitigs are arranged by starting position, not ending
  for ( i = 0; i < numUnitigs; i++)
  {
  	IntUnitigPos *upos = GetIntUnitigPos( ma->u_list, i);
	unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
	// int isSurrogate = unitig->flags.bits.isSurrogate; // (unitig->info.CI.baseID > 0);
	
	if (debug == 1)
	  fprintf( stderr, "in contig %d, unitig %d is at %f -> %f maxContigPos: %f, isSurrogate: %d, baseID: %d\n", 
			   contig->id, unitig->id,
			   unitig->offsetAEnd.mean, unitig->offsetBEnd.mean, maxContigPos,
			   isSurrogate, unitig->info.CI.baseID );

	if ( unitig->offsetAEnd.mean >= maxContigPos || unitig->offsetBEnd.mean >= maxContigPos )
	{
	  maxContigPos = max( unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
	  *unitigID = unitig->id;
	  isSurrogate = unitig->flags.bits.isSurrogate;
	}
  }
  
  if ( isSurrogate )
  {
	fprintf( stderr, "unitig %d on high end of contig %d is surrogate!\n",
			 unitig->id, contig->id);
	return 1;
  }
  else
	return 0;
}

// findFirstUnitig looks for a surrogate at the low end of a contig
int findFirstUnitig( ContigT *contig, int *unitigID )
{
  MultiAlignT *ma;
  int numUnitigs;
  IntUnitigPos *upos;
  NodeCGW_T *unitig;
  int isSurrogate;
  
  fprintf( stderr, "in FindFirstUnitig\n");
  
  ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numUnitigs = GetNumIntUnitigPoss( ma->u_list );
  
  upos = GetIntUnitigPos( ma->u_list, 0);
  unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
  isSurrogate = unitig->flags.bits.isSurrogate; // (unitig->info.CI.baseID > 0);
  
  if (debug == 1)
	  fprintf( stderr, "in contig %d, unitig %d is at %f -> %f, isSurrogate: %d\n", 
			   contig->id, unitig->id,
			   unitig->offsetAEnd.mean, unitig->offsetBEnd.mean,
			   isSurrogate);
  
  *unitigID = unitig->id;
  if ( isSurrogate )
	  // if ( unitig->offsetAEnd.mean == 0.0 || unitig->offsetBEnd.mean == 0.0 )
  {
	  fprintf( stderr, "unitig %d on low end of contig %d is surrogate!\n",
			   unitig->id, contig->id);
	  return 1;
  }
  else
	  return 0;
}

int examineGap( ContigT *lcontig, int lFragIid, ContigT *rcontig, int rFragIid, 
				int gapNumber, int *ahang, int *olapLengthOut, int *bhang, int *currDiffs,
				int *lcontigBasesIntact, int *rcontigBasesIntact,
				int *closedGapDelta, int lBasesToNextFrag, int rBasesToNextFrag,
				int *leftFragFlapLength, int *rightFragFlapLength)
{
  CIFragT *lFrag = NULL, *rFrag = NULL;
  static VA_TYPE(char) *ungappedSequence=NULL, *ungappedQuality=NULL;
  char lFragSeqBuffer[AS_BACTIG_MAX_LEN+1], lqltbuffer[AS_BACTIG_MAX_LEN+1];
  char rFragSeqBuffer[AS_BACTIG_MAX_LEN+1], rqltbuffer[AS_BACTIG_MAX_LEN+1];
  char lcompBuffer[AS_BACTIG_MAX_LEN+1], rcompBuffer[AS_BACTIG_MAX_LEN+1];
  uint lclr_bgn, lclr_end;
  uint rclr_bgn, rclr_end;
  char tmp_char, *lSequence, *rSequence;
  NodeOrient lContigOrientation, rContigOrientation;
  InfoByIID *info;
  int temp, len;
  int i;
  int lcontigBaseStart, lcontigBasesUsed, rcontigBasesUsed;
  int lFragContigOverlapLength, rFragContigOverlapLength;
  
  if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
	lContigOrientation = A_B;
  else
	lContigOrientation = B_A;
  
  if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
	rContigOrientation = A_B;
  else
	rContigOrientation = B_A;

  if ( lFragIid != -1)
  {
	info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, lFragIid);
	assert( info->set );
	lFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
  }
  
  if ( rFragIid != -1)
  {
	info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, rFragIid);
	assert( info->set );
	rFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
  }
  
  if( fsread == NULL )
  {
    fsread = new_ReadStruct();
  }
  if (ungappedSequence== NULL ) 
  {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality = CreateVA_char(0);
  } 
  else 
  {
    ResetVA_char(ungappedSequence);
    ResetVA_char(ungappedQuality);
  }

  if ( lFragIid != -1)
  {
	getFragStore( ScaffoldGraph->fragStore, lFragIid, FRAG_S_ALL, fsread);
	getClearRegion_ReadStruct( fsread, &lclr_bgn, &lclr_end, READSTRUCT_CNS);
	getSequence_ReadStruct( fsread, lFragSeqBuffer, lqltbuffer, AS_BACTIG_MAX_LEN);
  }
  
  if ( rFragIid != -1)
  {
	getFragStore( ScaffoldGraph->fragStore, rFragIid, FRAG_S_ALL, fsread);
	getClearRegion_ReadStruct( fsread, &rclr_bgn, &rclr_end, READSTRUCT_CNS);
	getSequence_ReadStruct( fsread, rFragSeqBuffer, rqltbuffer, AS_BACTIG_MAX_LEN);
  }
  
  // Get the consensus sequences for both chunks from the Store
  GetConsensus( ScaffoldGraph->ContigGraph, lcontig->id, lContigConsensus, lContigQuality);
  GetConsensus( ScaffoldGraph->ContigGraph, rcontig->id, rContigConsensus, rContigQuality);
  
  lSequence = Getchar( lContigConsensus, 0);
  rSequence = Getchar( rContigConsensus, 0);
  
  // ----------------------> lContigOrientation == A_B
  //                  -----> frag is 5p->3p into gap, aligned with contig
  
  // <---------------------- lContigOrientation == B_A
  //                  -----> frag is 5p->3p into gap, aligned opposite to contig
  
  if (lContigOrientation == B_A)  // the frag is oriented opposite to the contig in this case
  {
	SequenceComplement( lSequence, NULL);
  }
  
  // print out info on the left contig and fragment
  {
	if ( lFragIid != -1)
	{
	  char tmp_char = lFragSeqBuffer[lclr_end];
	  
	  lFragSeqBuffer[lclr_end] = '\0';
	  fprintf( stderr, " last 50 bases of lfragIid clr range: %s\n", &lFragSeqBuffer[lclr_end - 50]);
	  lFragSeqBuffer[lclr_end] = tmp_char;
	  
	  fprintf( stderr, "for frag %d, lclr_bgn: %d, lclr_end: %d, strlen( lFragSeqBuffer ): " F_SIZE_T "\n", 
			   lFragIid, lclr_bgn, lclr_end, strlen( lFragSeqBuffer ));
	  // fprintf( stderr, " lfrag: %s\n", lFragSeqBuffer);
	}
	
	fprintf( stderr, "  last 50 bases of lContig consensus: %s\n", &lSequence[ strlen( lSequence ) - 50]);
	fprintf( stderr, "\n");
  }
  
  // ----------------------> rContigOrientation == A_B
  // <-----                  frag is 5p->3p into gap, aligned opposite to contig
  
  // <---------------------- rContigOrientation == B_A
  // <-----                  frag is 5p->3p into gap, aligned with contig
  
  // now do right contig
  if (rContigOrientation == B_A)  // the frag is oriented opposite to the contig in this case
  {
	SequenceComplement( rSequence, NULL);  // flip contig sequence to its orientation in scaffold
  }
  
  if ( rFragIid != -1)
  {
	// we want to flip the frag in either case
	SequenceComplement( rFragSeqBuffer, NULL);
	len = strlen( rFragSeqBuffer );
	temp = len - rclr_bgn;  // new rclr_end
	rclr_bgn = len - rclr_end;
	rclr_end = temp;
  }
  
  // print out info on the right contig and fragment
  {
	if ( rFragIid != -1)
	{
	  char tmp_char = rFragSeqBuffer[rclr_bgn + 50];
	  rFragSeqBuffer[rclr_bgn + 50] = '\0';
	  fprintf( stderr, "first 50 bases of rFragIid clr range: %s\n", &rFragSeqBuffer[rclr_bgn]);
	  rFragSeqBuffer[rclr_bgn + 50] = tmp_char;
	  
	  fprintf( stderr, "for frag %d, rclr_bgn: %d, rclr_end: %d, strlen( rFragSeqBuffer ): " F_SIZE_T "\n", 
			   rFragIid, rclr_bgn, rclr_end, strlen( rFragSeqBuffer ));
	  // fprintf( stderr, " rfrag: %s\n", rFragSeqBuffer);
	}
	
	tmp_char = rSequence[50];
	rSequence[50] = '\0';
	fprintf( stderr, " first 50 bases of rContig consensus: %50s\n", rSequence);
	rSequence[50] = tmp_char;
	
	fprintf( stderr, "\n");
  }
  
  // we use frag sequence from where the clear range ends to the end of the frag
  // ----------------------> lContigOrientation == A_B
  //               ----->    frag is 5p->3p into gap, aligned with contig
  // <---------------------- lContigOrientation == B_A
  //               ----->    frag is 5p->3p into gap, aligned opposite to contig
  
  if ( lFragIid != -1)
  {
	if ( lContigOrientation == A_B )
	  lFragContigOverlapLength = (int) (lcontig->bpLength.mean - lFrag->contigOffset3p.mean);
	else
	  lFragContigOverlapLength = (int) (lFrag->contigOffset3p.mean);
  }
  else
	lFragContigOverlapLength = 0;

  // grab the last CONTIG_BASES bases of the lcontig consensus sequence
  lcontigBasesUsed = min( CONTIG_BASES - lFragContigOverlapLength, 
						  (int) lcontig->bpLength.mean - lFragContigOverlapLength);
  // lcontigBasesUsed = 100.0;  // temp hack, but it is sometimes better to do 100 than 1000.  why???
  
  // lcontigBaseStart is the base where we start using the consensus sequence in lcompBuffer
  // and thus also the number of bases from the contig that are intact
  lcontigBaseStart = strlen( lSequence ) - lcontigBasesUsed - lFragContigOverlapLength;
  
  // grab the bases from the contig, ie, those not from the non-clear range of the frag
  for ( i = 0; i < lcontigBasesUsed; i++)
  {
	lcompBuffer[ i ] = lSequence[ lcontigBaseStart + i];  // a bit ugly
  }
  lcompBuffer[ i ] = '\0';

  // now tack on the 3p clr range extension to the bases of the contig consensus sequence
  if ( lFragIid != -1) // && 0)  // temp hack
  {
	// if (lcontigBasesUsed < lFragContigOverlapLength)  // this means that the frag and contig overlap by
	// lFragContigOverlapLength = lcontigBasesUsed;    // more bases than we are using from the contig


	// basesToNextFrag is the number of bases back to the first frag that gets us to 2x
      // used to be, but Aaron thought it looked fishy
	//MaxEndGap = strlen( lFragSeqBuffer ) - lclr_end - lFragContigOverlapLength + lBasesToNextFrag + 20;  // 20 is slop
	MaxEndGap = strlen( lFragSeqBuffer ) - lclr_end + lBasesToNextFrag + 20;  // 20 is slop
	fprintf(stderr,"## MaxEndGap %d\n",MaxEndGap);
	for ( i = lclr_end; i < strlen( lFragSeqBuffer ); i++)
	  lcompBuffer[ lcontigBasesUsed + i - lclr_end ] = lFragSeqBuffer[ i ];
	lcompBuffer[ lcontigBasesUsed + i - lclr_end ] = '\0';
  }
  else
  {
	MaxEndGap = 100;
  }
  
  // we use frag sequence from where the clear range ends to the end of the frag
  // ----------------------> rContigOrientation == A_B
  //    <-----               frag is 5p->3p into gap, aligned opposite to contig

  // <---------------------- rContigOrientation == B_A
  //      <-----             frag is 5p->3p into gap, aligned with contig

  if (rFragIid != -1)
  {
	if ( rContigOrientation == A_B )
	  rFragContigOverlapLength = (int) (rFrag->contigOffset3p.mean);
	else
	  rFragContigOverlapLength = (int) (rcontig->bpLength.mean - rFrag->contigOffset3p.mean);
  }
  else
	rFragContigOverlapLength = 0;
  
  if ( rFragIid != -1)
  {
	// basesToNextFrag is the number of bases back to the first frag that gets us to 2x
      // used to be, but Aaron thought it looked fishy
	//MaxBegGap = rclr_bgn - rFragContigOverlapLength + rBasesToNextFrag + 20;  // 20 is slop	
	MaxBegGap = rclr_bgn + rBasesToNextFrag + 20;  // 20 is slop	
	fprintf(stderr,"## MaxBegGap %d\n",MaxBegGap);
  }
  else
  {
	MaxBegGap = 100;
  }
  
  // now if we have a right frag, grab the "5p" clr range extension - remember the frag has been flipped
  if ( rFragIid != -1) //  && 0)  // temp hack
  {
	for ( i = 0; i < rclr_bgn; i++)
	  rcompBuffer[ i ] = rFragSeqBuffer[ i ];
  }
  else
  {
	rclr_bgn = 0;  // need this for the next loop
  }
  
  // grab the first CONTIG_BASES bases of the rcontig consensus sequence
  // the rcontig consensus has been flipped if necessary
  rcontigBasesUsed = min( CONTIG_BASES - rFragContigOverlapLength, 
						  (int) rcontig->bpLength.mean - rFragContigOverlapLength);
  for ( i = 0; i < rcontigBasesUsed; i++)
	rcompBuffer[ rclr_bgn + i ] = rSequence[ i + rFragContigOverlapLength ]; //Aaron ad rFragContigOverlapLength
  rcompBuffer[ rclr_bgn + i ] = '\0';

  fprintf( stderr, "> lcompBuffer gap %d (len: " F_SIZE_T "): \n%s\n", gapNumber, strlen( lcompBuffer ), lcompBuffer);
  fprintf( stderr, "> rcompBuffer gap %d (len: " F_SIZE_T "): \n%s\n", gapNumber, strlen( rcompBuffer ), rcompBuffer);
  
  // now lcompBuffer and rcompBuffer hold the sequence of the fragments in the correct strand
  // now prepare for call to Local_Overlap_AS_forCNS
  {
	int beg, end, opposite = FALSE;
	double erate, thresh, minlen;
	CompareOptions what;
	Overlap *overlap;
	int basesAdded;
	LengthT gapSize;
	char *rcompBufferTrimmed = NULL;
	
	// MaxGaps = 5;  // we don't want a lot of segments
	
	// stole from MultiAlignment_CNS.c
#define CNS_DP_ERATE .12 // was .20 during testing
#define CNS_DP_THRESH 1e-6
#define CNS_DP_MINLEN 30

	beg = -strlen (rcompBuffer);
	end = strlen (lcompBuffer);
	erate = CNS_DP_ERATE;
	thresh = CNS_DP_THRESH;
	minlen = CNS_DP_MINLEN;
	// next line was: what = AS_FIND_LOCAL_ALIGN_NO_TRACE;
	what = AS_FIND_LOCAL_ALIGN;
	
	/* if ( ( min (rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean) - 
	   max (lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean) ) < -19.0) */ // part 1 of temp hack

	fprintf( stderr, "MaxBegGap: %d\n", MaxBegGap);
	fprintf( stderr, "MaxEndGap: %d\n", MaxEndGap);
	
	overlap = Local_Overlap_AS_forCNS( lcompBuffer, rcompBuffer,
									   -strlen( rcompBuffer) , strlen( lcompBuffer), opposite,
									   erate, thresh, minlen,
									   what);
	/* else
	   overlap = NULL; */ // part 2 of temp hack

	if (1 && overlap != NULL)
	{
	  fprintf( stderr, "initial ahang: %d, bhang:%d, length: %d, diffs: %d, diffs / length %%: %f\n",
			   overlap->begpos, overlap->endpos, overlap->length, overlap->diffs, 
			   100.0 * overlap->diffs / overlap->length);
	  
	  // printOverlap( overlap, lcompBuffer, rcompBuffer);  // my crappy routine
	  Print_Overlap( stderr, lcompBuffer, rcompBuffer, overlap);
	}
	
	// not interested in overlaps with negative ahangs or bhangs
	if (overlap != NULL && overlap->begpos > 0 && overlap->endpos > 0)
	{
          /*
	       \
	        \
	         \  left flap, right frag
	          \  
	           -------------------
	       -------------------
	                          \   right flap, left frag
	                           \
	                            \
          */
          
	  // do flap trimming
	  // left flap, right frag
	  *rightFragFlapLength = 0;
	  if ( overlap->begpos == abs(overlap->trace[0]) - 1 && overlap->trace[0] < 0 )
	  {
		while ( overlap->begpos == abs(overlap->trace[ *rightFragFlapLength ]) - 1 && 
				overlap->trace[ *rightFragFlapLength ] < 0 )
		  (*rightFragFlapLength)++;
	  }
	  
	  // right flap, left frag
	  {
	  	// first find the number of elements in trace
		int numTraceElements = 0, icnt = 0;
		while ( overlap->trace[ icnt++ ] != 0 ) numTraceElements++;
		
		*leftFragFlapLength = 0;
		if ( overlap->trace[ numTraceElements - 1 ] > 0 && overlap->endpos > 0)
		{
		  icnt = numTraceElements - 1;
		  while (overlap->trace[ icnt ] == strlen(rcompBuffer) - overlap->endpos + 1 && icnt >= 0 )
		  {
			(*leftFragFlapLength)++; icnt--;
		  }
		}
	  }
	
	  lcompBuffer[ strlen(lcompBuffer) - *leftFragFlapLength ] = '\0';
	  // MaxEndGap -= *leftFragFlapLength;
	  rcompBufferTrimmed = &rcompBuffer[ *rightFragFlapLength ];
	  // MaxBegGap = *rightFragFlapLength;
	  
	  // now do overlap again after trimming to make sure it is still there, sometimes trimming makes them go away
	  if (1)  // hack
		overlap = Local_Overlap_AS_forCNS( lcompBuffer, rcompBufferTrimmed,
										   -strlen( rcompBufferTrimmed) , strlen( lcompBuffer), opposite,
										   erate, thresh, minlen,
										   what);
	  if (!overlap)
		fprintf( stderr, "lost overlap to flap trimming!\n");
	}
	
	if (1 && overlap != NULL && overlap->begpos > 0 && overlap->endpos > 0)
	{
	  fprintf( stderr, "post-flap trimming ahang: %d, bhang:%d, length: %d, diffs: %d, diffs / length %%: %f\n",
			   overlap->begpos, overlap->endpos, overlap->length, overlap->diffs, 
			   100.0 * overlap->diffs / overlap->length);
	  
	  // printOverlap( overlap, lcompBuffer, rcompBuffer);  // my crappy routine
	  Print_Overlap( stderr, lcompBuffer, rcompBufferTrimmed, overlap);
	}
	
	// not interested in overlaps with negative ahangs or bhangs
	if (overlap != NULL && overlap->begpos > 0 && overlap->endpos > 0)
	{
	  int baseChangeLeftContig, baseChangeRightContig;

	  fprintf( stderr, "found overlap between frags %d and %d, length = %d\n", 
			   lFragIid, rFragIid, overlap->length);
	  fprintf( stderr, "ahang + overlap->length: %d, strlen( lcompBuffer ): " F_SIZE_T ", diff: " F_SIZE_T "\n",
			   overlap->begpos + overlap->length, strlen( lcompBuffer ),
			   overlap->begpos + overlap->length - strlen( lcompBuffer ));
	  fprintf( stderr, "overlap->length + bhang: %d, strlen( rcompBufferTrimmed ): " F_SIZE_T ", diff: " F_SIZE_T "\n",
			   overlap->length + overlap->endpos, strlen( rcompBufferTrimmed ), 
			   overlap->length + overlap->endpos - strlen( rcompBufferTrimmed ));
	  
	  *lcontigBasesIntact = lcontigBaseStart;
	  *ahang = max (0, overlap->begpos);
	  *olapLengthOut = overlap->length;
	  *bhang = max (0, overlap->endpos);
	  *currDiffs = overlap->diffs;
	  *rcontigBasesIntact = max( (int) rcontig->bpLength.mean - CONTIG_BASES, 0);

	  // calculate how many bases would be changed if gap was closed
	  // take the whole lcompBuffer, subtract lcontigBasesUsed and lFragContigOverlapLength from ahang
	  // and add in the length of the overlap
	  // take the whole rcompBuffer, subtract rcontigBasesUsed and rFragContigOverlapLength from bhang

	  baseChangeLeftContig = overlap->begpos - lcontigBasesUsed - lFragContigOverlapLength;
	  baseChangeRightContig = overlap->endpos - rcontigBasesUsed - rFragContigOverlapLength;

	  fprintf( stderr, "lcontigBasesIntact: %d\n", *lcontigBasesIntact);
	  fprintf( stderr, "overlap->begpos: %d\n", overlap->begpos);
	  fprintf( stderr, "overlap->length: %d\n", overlap->length);
	  fprintf( stderr, "overlap->endpos: %d\n", overlap->endpos);
	  fprintf( stderr, "rcontigBasesIntact: %d\n", *rcontigBasesIntact);

	  // fprintf( stderr, "base change left contig: %d\n",
	  //   overlap->begpos - lcontigBasesUsed - lFragContigOverlapLength);
	  // fprintf( stderr, "base change right contig: %d\n",
	  //   overlap->endpos - rcontigBasesUsed - rFragContigOverlapLength);	  
	  
	  basesAdded = baseChangeLeftContig + overlap->length + baseChangeRightContig;
	  
	  gapSize = FindGapLength( lcontig, rcontig, FALSE);

	  totalContigsBaseChange += basesAdded;
	  
	  *closedGapDelta = basesAdded - (int) gapSize.mean;	  
	  
	  fprintf( stderr, "would fill gap %d of size %d with %d bases, net change: %d\n",
			   gapNumber, (int) gapSize.mean, basesAdded, basesAdded - (int) gapSize.mean);
	  fprintf( stderr, "lcontig->bpLength.mean: %f, baseChangeLeftContig: %d\n",
			   lcontig->bpLength.mean, baseChangeLeftContig);
	  fprintf( stderr, "rcontig->bpLength.mean: %f, baseChangeRightContig: %d\n",
			   rcontig->bpLength.mean, baseChangeRightContig);

	  fprintf( stderr, "new contig size: %d\n", 
			   (int) lcontig->bpLength.mean + baseChangeLeftContig + 
			   (int) rcontig->bpLength.mean + baseChangeRightContig + overlap->length);
	  fprintf( stderr, "totalContigsBaseChange: %d\n", totalContigsBaseChange);
	  
	  // *closedGapSizeDiff = basesAdded - (int) gapSize.mean;

	  return 1;
	}
	else
	{
	  fprintf( stderr, "no overlap found between frags %d and %d\n", lFragIid, rFragIid);
	  return 0;
	}
  }
}

int dumpFragInfo( int fragIid )
{
  CIFragT *frag;
  InfoByIID *info;    

  info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, fragIid);
  assert(info->set);
  frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
  
  fprintf( stderr, "    frag: %6d\t contig id: %d\t contig pos (5p, 3p): %6d, %6d\n", 
	       fragIid, frag->contigID, (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);	  
  
  return 0;
}

void printOverlap( Overlap *olap, char *aseq, char *bseq )
{
  int i, j;
  int ahang = olap->begpos;
  int bhang = olap->endpos;
  int length = olap->length;
  int numLines = 1 + (ahang + length + bhang) / 50;
  int astart, bstart, aend, bend;

#if 0
  { 
	int *trace;
	OverlapMesg *align, alignee;

	align = &alignee;
	
	align->ahg = olap->begpos;
	align->delta = olap->trace;

	trace = my_Unpack_Alignment_AS(align);
	
	fprintf( stderr, "\nUncompressed trace:\n");
	for (i = 0; trace[i] != 0; i++)
	  fprintf( stderr, "myTrace:  %3d\n", trace[i]);

	i = 0;
	while (olap->trace[i] != 0)
	{
	  fprintf( stderr, "trace[ %d ] = %d\n", i, olap->trace[i]);
	  i++;
	}
  }
#endif
  
  if (ahang > 0)
  {
	astart = 0;
	bstart = ahang;
	aend = ahang + length;
	bend = bstart + length;
  }
  else
  {
	astart = -ahang;
	bstart = 0;
	aend = astart + length;
	bend = -ahang + bstart + length;
  }
  if (bhang > 0)
	bend += bhang;
  else
	aend -= bhang;
  

  fprintf( stderr, "astart: %d, aend: %d\n", astart, aend);
  fprintf( stderr, "bstart: %d, bend: %d\n", bstart, bend);  
  
  for ( i = 0; i < numLines; i++)
  {
	// do aseq
	for (j = 0; j < 50; j++)
	{
	  if (i * 50 + j >= astart && i * 50 + j < aend)
		fprintf( stderr, "%c", aseq[ i * 50 + j - astart]);
	  else
		fprintf( stderr, " ");
	}
	fprintf( stderr, " lcomp\n");
	
	// do bseq
	for (j = 0; j < 50; j++)
	{
	  if (i * 50 + j >= bstart && i * 50 + j < bend)
		fprintf( stderr, "%c", bseq[ i * 50 + j - bstart]);
	  else
		fprintf( stderr, " ");
	}
	fprintf( stderr, " rcomp\n");
	fprintf( stderr, "\n");
  }
}

#if 0
	// diagnostic info, stolen from AS_ALN_testdriver5.c
	if (1 && overlap != NULL)
	{
	  InternalFragMesg A, B;
	  int alen, blen, del,sub, ins, affdel, affins, blockdel, blockins;
	  double errRate, errRateAffine;
	  int AFFINEBLOCKSIZE=4;
	  
	  A.sequence = lcompBuffer;
	  B.sequence = rcompBuffer;
	  
	  fprintf( stderr, "about to Print_Overlap_AS\n");
	  fflush ( stderr );
	  
	  Print_Overlap_AS( stdout, &A, &B, overlap);
	  fflush ( stdout );
	  
	  Analyze_Affine_Overlap_AS(&A, &B, overlap, AS_ANALYZE_ALL, &alen, &blen, &del, &sub, &ins,
								&affdel, &affins, &blockdel, &blockins, AFFINEBLOCKSIZE);
	  
	  errRate = (sub + ins + del) / (double)(alen + ins);
	  
	  errRateAffine = (sub + affins + affdel) /
		(double)(alen - del + affins + affdel);
	  
	  printf("\n\nAlen %d, Blen %d, del %d, sub %d, ins %d\n"
			 " affdel %d, affins %d, blockdel %d, blockins %d\n",
			 alen, blen, del, sub, ins,
			 affdel, affins, blockdel, blockins);
	  printf("Simple mismatch rate %f\n", errRate);
	  printf("Affine mismatch rate %f\n", errRateAffine);
	  fflush ( stdout );
	}
#endif

int *my_Unpack_Alignment_AS(OverlapMesg *align)
{ static int    UnpackBuffer[2*AS_READ_MAX_LEN+1];
  signed char  *calign;
  int           apos, bpos;

  apos   = align->ahg + 1;  /* ahg >= 0 for all overlaps */
  bpos   = 1;
  calign = align->delta;

  { int i, c;
    int *spt;

    spt = UnpackBuffer;
    for (i = 0; (c = calign[i]) != 0; i++)
      if (c == AS_LONG_DELTA_CODE)
        { apos += AS_LONGEST_DELTA;      /* Uninterrupted match of 126 bases */
          bpos += AS_LONGEST_DELTA;
        }
      else if (c == AS_POLY_DELTA_CODE)
        { c = calign[++i];  /* Multi-base gap */
          if (c < 0)
            while (c < 0)
              { c    += 1;
                bpos += 1;
                *spt++ = -apos;
              }
          else
            while (c > 0)
              { c    -= 1;
                apos += 1;
                *spt++ = -bpos;
              }
        }
      else
        { if (c < 0)        /* Single gap */
            { bpos -= c;
              apos -= (c+1);
              *spt++ = -apos;
            }
          else
            { bpos += (c-1);
              apos += c;
              *spt++ = bpos;
            }
        }
    *spt = 0;
  }

  return (UnpackBuffer);
}

void extendContig( ContigT *contig, int extendAEnd)
{
  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd
  int contigOrientation = contig->offsetAEnd.mean < contig->offsetBEnd.mean ? A_B : B_A;
  MultiAlignT *cma = CreateEmptyMultiAlignT();
  float lengthDelta;
  
  cma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE);

  lengthDelta = strlen( Getchar( cma->consensus, 0 )) - contig->bpLength.mean;

  if (contigOrientation == A_B)
	fprintf( stderr, "in extendContig, contig %8d original pos: %.0f, %.0f\n",
			 contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
  else
	fprintf( stderr, "in extendContig, contig %8d original pos: %.0f, %.0f\n",
			 contig->id, contig->offsetBEnd.mean, contig->offsetAEnd.mean);
  
  contig->bpLength.mean += lengthDelta;
  
  if ( contigOrientation == A_B )
  {
	if ( extendAEnd == TRUE )
	  contig->offsetAEnd.mean -= lengthDelta;
	else  // extend the B end
	  contig->offsetBEnd.mean += lengthDelta;
  }
  else  // contigOrientation is B_A
  {
	if ( extendAEnd == TRUE )
	  contig->offsetAEnd.mean += lengthDelta;
	else  // extend the B end
	  contig->offsetBEnd.mean -= lengthDelta;
  }

  if (contigOrientation == A_B)
	fprintf( stderr, "in extendContig, contig %8d altered pos: %.0f, %.0f\n",
			 contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
  else
	fprintf( stderr, "in extendContig, contig %8d altered pos: %.0f, %.0f\n",
			 contig->id, contig->offsetBEnd.mean, contig->offsetAEnd.mean);  
}

void extendUnitigs( NodeCGW_T *unitig, int fragIid, extendableFrag extFrag, int leftContig)
{
  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd
  int unitigLeftEnd, unitigRightEnd, unitigOrientation;
  int AEnd;
  
  // GetContigPositionInScaffold( contig, &contigLeftEnd, &contigRightEnd, &contigOrientation);
  GetChunkPositionInContig( unitig, &unitigLeftEnd, &unitigRightEnd, &unitigOrientation);

  unitig->bpLength.mean += extFrag.extension;
  
  // the unitig we're extending has to be on one end or the other of the contig
  if ( unitigLeftEnd == 0 || unitigRightEnd == 0)
	AEnd = TRUE;
  else
	AEnd = FALSE;

  if ( AEnd == TRUE )
  {
	ContigTIterator contigIter;
	NodeCGW_T *currUnitig;
	ContigT *contig = GetGraphNode( ScaffoldGraph->ContigGraph, unitig->info.CI.contigID);
	
	InitContigTIterator( ScaffoldGraph, contig->id, TRUE, FALSE, &contigIter);
	while ( NextContigTIterator( &contigIter ))
	{
	  currUnitig = GetGraphNode( ScaffoldGraph->CIGraph, contigIter.curr);

	  if (currUnitig->id != unitig->id)  // adjust everybody upwards except extended unitig
	  {
		currUnitig->offsetAEnd.mean += extFrag.extension;
		currUnitig->offsetBEnd.mean += extFrag.extension;
	  }
	  else
	  {
		if ( unitigOrientation == 0 /* A_B */)
		  currUnitig->offsetBEnd.mean += extFrag.extension;
		else  // unitig is B_A in contig
		  unitig->offsetAEnd.mean += extFrag.extension;
	  }
	}
  }
  else
  {
	if ( unitigOrientation == 0 /* A_B */)
	{
	  fprintf (stderr, "extending AEnd of unitig %8d from %.0f to %.0f\n",
			   unitig->id, unitig->offsetAEnd.mean, unitig->offsetAEnd.mean + extFrag.extension);
	  unitig->offsetBEnd.mean += extFrag.extension;
	}
	else // unitigOrientation == 1 /* B_A */
	{
	  fprintf (stderr, "extending AEnd of unitig %8d from %.0f to %.0f\n",
			   unitig->id, unitig->offsetAEnd.mean, unitig->offsetAEnd.mean + extFrag.extension);
	  unitig->offsetAEnd.mean += extFrag.extension;
	}
  }
}

// stole most of this from OutputUnitigsFromMultiAligns
/********************************************************************************/
int GetNewUnitigMultiAlign( NodeCGW_T *unitig, fragPositions *fragPoss, int extendedFragIid)
{
  GenericMesg			pmesg;
  IntUnitigMesg			ium_mesg;
  int i;
  int numCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  UnitigStatus   status;
  // IUMStruct * is;
  
#define  USE_UNGAPPED_CONSENSUS_FOR_UNITIG
#ifdef USE_UNGAPPED_CONSENSUS_FOR_UNITIG
  static VA_TYPE(char) *ungappedSequence=NULL,*ungappedQuality=NULL;
  if (ungappedSequence== NULL ) {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality = CreateVA_char(0);
  } else {
    ResetVA_char(ungappedSequence);
    ResetVA_char(ungappedQuality);
  }
#endif


  pmesg.m = &ium_mesg;
  pmesg.t = MESG_IUM;

  //    MultiAlignT *ma = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, unitig->id);

  assert( unitig->id>=0 && unitig->id< numCIs);

  switch(unitig->type)
  {
    case DISCRIMINATORUNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      //	fprintf(GlobalData->stderrc,"* Unitig %d is UNIQUE: DISCRIMINATOR output %d \n",unitig->id, cid);
      break;
    case UNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      //	fprintf(GlobalData->stderrc,"* Unitig %d is UNIQUE: output %d \n",unitig->id, cid);
      break;
    case UNRESOLVEDCHUNK_CGW:
      if(unitig->info.CI.numInstances > 0)
	  {
		assert(!unitig->flags.bits.isUnique);
		status = AS_SEP;
		//	fprintf(GlobalData->stderrc,
		//  "* Unitig %d has %d instances--- output %d SEP\n",unitig->id, unitig->info.CI.numInstances,cid);
      }
	  else
	  {
		if(unitig->scaffoldID != NULLINDEX)
		{
		  //	  fprintf(GlobalData->stderrc, 
		  //      "* Unitig %d has %d instances--- output %d UNIQUE\n",unitig->id, unitig->info.CI.numInstances,cid);
		  status = AS_UNIQUE;
		}
		else
		{
		  //	  fprintf(GlobalData->stderrc,
		  //      "* Unitig %d has %d instances--- output %d NOTREZ\n",unitig->id, unitig->info.CI.numInstances,cid);
		  status = AS_NOTREZ;
		}
      }
      break;
    case RESOLVEDREPEATCHUNK_CGW:
      /* SKIP THESE */
      //      fprintf(GlobalData->stderrc,"* Skipping unitig %d --- RESOLVEDREPEAT\n",unitig->id);
      // continue;
    default:
      assert(0);
  }

  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, unitig->id, TRUE);

  //fprintf( stderr, "for unitig %d, before reforming, strlen( ma->consensus ) = " F_SIZE_T "\n",
  //	   unitig->id, strlen( Getchar( ma->consensus, 0) ));
  //fprintf( stderr, "for unitig %d, before reforming, consensus:\n %s\n",
  //	   unitig->id, Getchar( ma->consensus, 0) );

  {
    // int numFrag = GetNumIntMultiPoss(ma->f_list);
    int extendedFragLeftward, aligned;
    // int32 tmpSource[GetNumIntMultiPoss(ma->f_list) + 1];
    assert (unitig->type != CONTIG_CGW);

#if 0
    {
      int i;
      // Null out the source field
      for( i = 0; i < numFrag; i++)
	  {
		IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
		tmpSource[i] = (int) mp_i->source;
		mp_i->source = NULL;
		assert(mp_i->ident);
      }
    }
#endif
    
    ium_mesg.iaccession = unitig->id;

    ium_mesg.source = NULL;

    ium_mesg.coverage_stat = unitig->info.CI.coverageStat;
    ium_mesg.status = status;
    ium_mesg.a_branch_point = unitig->info.CI.branchPointA;
    ium_mesg.b_branch_point = unitig->info.CI.branchPointB;

#ifdef USE_UNGAPPED_CONSENSUS_FOR_UNITIG
    GetMultiAlignUngappedConsensus(ma, ungappedSequence, ungappedQuality);
    ium_mesg.consensus = Getchar(ungappedSequence,0);
    ium_mesg.quality = Getchar(ungappedQuality,0);
    ium_mesg.length = GetMultiAlignUngappedLength(ma);
#else
    ium_mesg.length = GetMultiAlignLength(ma);
    ium_mesg.consensus = Getchar(ma->consensus,0);
    ium_mesg.quality = Getchar(ma->quality,0);
#endif
    ium_mesg.forced = 0;
    ium_mesg.num_frags = GetNumIntMultiPoss(ma->f_list);

	// replace the positions in the f_list with the adjusted positions
	extendedFragLeftward = FALSE;
	for ( i = 0; i < GetNumIntMultiPoss( ma->f_list ); i++)
	{
	  IntMultiPos *tempPos = GetIntMultiPos( ma->f_list, i);

	  // fprintf( stderr, "in GetNewUnitigMultiAlign, changing frag %10d from ( %8d, %8d ) to ( %8d, %8d )\n",
	  //   tempPos->ident, tempPos->position.bgn, tempPos->position.end, fragPoss[i].bgn, fragPoss[i].end);

	  if (tempPos->ident == extendedFragIid)  
	  {
		tempPos->contained = 0;    // the extended frag is no longer contained by any other frag
		if ( fragPoss[i].bgn == 0 || fragPoss[i].end == 0)
		  extendedFragLeftward = TRUE;  // if at the beginning of the unitig have to reorder frags
	  }
	  tempPos->position.bgn = fragPoss[i].bgn;
	  tempPos->position.end = fragPoss[i].end;
	}
    ium_mesg.f_list = GetIntMultiPos( ma->f_list, 0);
	
	if ( extendedFragLeftward )  // definitely reorder frags in f_list if we extended a frag leftward
	  leftShiftIUM( ium_mesg.f_list, GetNumIntMultiPoss( ma->f_list ), extendedFragIid);
	else  // might need to reorder frags in f_list if we extended a frag rightward
	  rightShiftIUM( ium_mesg.f_list, GetNumIntMultiPoss( ma->f_list ), extendedFragIid);

    // (GlobalData->writer)(GlobalData->outfp,&pmesg);
#if 0
    {
      int i;
      // Restore the source field
      for(i = 0; i < numFrag; i++){
		IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
		mp_i->source = (char *)tmpSource[i];
      }
    }
#endif

	aligned = MultiAlignUnitig( &ium_mesg, 
								ScaffoldGraph->fragStore,
								reformed_consensus,
								reformed_quality,
								reformed_deltas,
								CNS_STATS_ONLY,
								1,
								Local_Overlap_AS_forCNS);   // DP_Compare);
	if ( aligned == -1 )
	{
	  fprintf( stderr, "MultiAlignUnitig failure on unitig %d\n", unitig->id);
	  return FALSE;   // assert(0);
	}
  
	fprintf( stderr, "for unitig %d, after reforming, strlen( reformed_consensus ) = " F_SIZE_T "\n",
			 unitig->id, strlen( Getchar( reformed_consensus, 0) ));

	{
	  // MultiAlignT *new_ma = CreateMultiAlignTFromIUM( &ium_mesg, GetNumCIFragTs( ScaffoldGraph->CIFrags ), FALSE);
	  // "-2" tells CreateMultiAlignTFromIUM to copy over the source field in an IntMultiPos from the ium_mesg
	  MultiAlignT *new_ma = CreateMultiAlignTFromIUM( &ium_mesg, -2, FALSE);

	  // if (new_ma == NULL) ...   handle failure, but must always unload multialign since we changed it

	  // This adds a reference only if keepInCache is true...
	  UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE);
	  InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE, new_ma, FALSE);

	  // PrintMultiAlignT( GlobalData->stderrc, new_ma, ScaffoldGraph->fragStore, 0x0, 0, TRUE, 0, READSTRUCT_LATEST);

	  // DuplicateEntryInSequenceDB( ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE, ium_mesg.iaccession, TRUE, TRUE);
	  // UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);
	}
  }

  // DeleteMultiAlignT(ma);
  fflush(NULL);
  return TRUE;
}

int extendCgwClearRange( int fragIid, int frag3pDelta)
{
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;

  if (frag3pDelta < 0)
	fprintf( stderr, "Warning: frag3pDelta less than zero: %d\n", frag3pDelta);
  
  if( fsread == NULL ) fsread = new_ReadStruct();

#ifdef OLD_WAY
  if ( fragIid != -1)
  {
	FragStoreHandle fragStore;
	char buffer[256];
    sprintf(buffer,"%s", GlobalData->Frag_Store_Name);

	// set and close the current frag store
	closeFragStore( ScaffoldGraph->fragStore );
	
	// open it read/write and set the fragment
	ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw+");
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	setClearRegion_ReadStruct( fsread, clr_bgn, clr_end + frag3pDelta, READSTRUCT_CGW);
	setStatus = setFragStore( ScaffoldGraph->fragStore, fragIid, fsread);
	closeFragStore( ScaffoldGraph->fragStore );
	
	// open it as before
	ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw");
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	

	fprintf( stderr, "extendCgwClearRange, changed frag %d clr_end from %d to %d\n",
			 fragIid, clr_end, clr_end + frag3pDelta);
  }

#else

  if ( fragIid != -1)
  {
	// FragStoreHandle fragStore;
	// char buffer[256];
    // sprintf(buffer,"%s", GlobalData->Frag_Store_Name);

	// set and close the current frag store
	// closeFragStore( ScaffoldGraph->fragStore );
	
	// open it read/write and set the fragment
	// ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw+");
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	setClearRegion_ReadStruct( fsread, clr_bgn, clr_end + frag3pDelta, READSTRUCT_CGW);
	setStatus = setFragStore( ScaffoldGraph->fragStore, fragIid, fsread);
	// closeFragStore( ScaffoldGraph->fragStore );
	
	// open it as before
	// ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw");
	// getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	

	fprintf( stderr, "extendCgwClearRange, changed frag %d clr_end from %d to %d\n",
			 fragIid, clr_end, clr_end + frag3pDelta);
  }
#endif

  return ( setStatus ); 
}

int revertToCnsClearRange_old( int fragIid )
{
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	char buffer[256];
    sprintf(buffer,"%s", GlobalData->Frag_Store_Name);

	// set and close the current frag store
	closeFragStore( ScaffoldGraph->fragStore );
	
	// open it read/write and set the fragment
	ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw+");
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	setClearRegion_ReadStruct( fsread, clr_bgn, clr_end, READSTRUCT_CGW);
	setStatus = setFragStore( ScaffoldGraph->fragStore, fragIid, fsread);
	closeFragStore( ScaffoldGraph->fragStore );
	
	// open it as before
	ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw");
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);
  }
  return ( setStatus ); 
}

int revertToCnsClearRange( int fragIid )
{
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	// FragStoreHandle fragStore;
	// char buffer[256];
    // sprintf(buffer,"%s", GlobalData->Frag_Store_Name);

	// set and close the current frag store
	// closeFragStore( ScaffoldGraph->fragStore );
	
	// open it read/write and set the fragment
	// ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw+");
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	setClearRegion_ReadStruct( fsread, clr_bgn, clr_end, READSTRUCT_CGW);
	setStatus = setFragStore( ScaffoldGraph->fragStore, fragIid, fsread);
	// closeFragStore( ScaffoldGraph->fragStore );
	
	// open it as before
	// ScaffoldGraph->fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw");
	// getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);
  }
  return ( setStatus ); 
}

int printClearRanges( int fragIid )
{
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	FragStoreHandle fragStore;

	fragStore = openFragStore( GlobalData->Frag_Store_Name, "rw");
	getFragStore( fragStore, fragIid, FRAG_S_ALL, fsread);	

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_ORIGINAL);
	fprintf( stderr, "printClearRanges, frag %8d READSTRUCT_ORG clr_bgn: %5d, clr_end %5d, len: %4d\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_OVL);
	fprintf( stderr, "printClearRanges, frag %8d READSTRUCT_OVL clr_bgn: %5d, clr_end %5d, len: %4d\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	fprintf( stderr, "printClearRanges, frag %8d READSTRUCT_CNS clr_bgn: %5d, clr_end %5d, len: %4d\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CGW);
	fprintf( stderr, "printClearRanges, frag %8d READSTRUCT_CGW clr_bgn: %5d, clr_end %5d, len: %4d\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_LATEST);
	fprintf( stderr, "printClearRanges, frag %8d READSTRUCT_LAT clr_bgn: %5d, clr_end %5d, len: %4d\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	closeFragStore( fragStore );
  }

  return ( setStatus ); 
}

// are things in the frag store always stored 5p -> 3p?
int setCgwClearRange( int fragIid, int frag3pDelta)
{
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	setClearRegion_ReadStruct( fsread, clr_bgn, clr_end + frag3pDelta, READSTRUCT_CGW);
	
	setStatus = setFragStore( ScaffoldGraph->fragStore, fragIid, fsread);

	fprintf( stderr, "setCgwClearRange, changed frag %d clr_end from %d to %d\n",
			 fragIid, clr_end, clr_end + frag3pDelta);
  }
  return ( setStatus ); 
}

int setCgwClearRangeAbs( int fragIid, int new_clr_bgn, int new_clr_end)
{
  unsigned int clr_bgn, clr_end;
  int offset, setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);

	fprintf( stderr, "as read, frag %d clr_bgn %d, clr_end %d\n",
			 fragIid, clr_bgn, clr_end);

	setClearRegion_ReadStruct( fsread, new_clr_bgn, new_clr_end, READSTRUCT_CGW);
	
	{
	  FILE *srcStore = fopen( "chr21_newer.frgStore/db.frg", "r+");
	  
	  if ( srcStore == NULL)
	  {
		fprintf( stderr, "failed to open srcStore\n");
		exit(1);
	  }
	  offset = fragIid * 56 - 8;
	  fprintf( stderr, "offset: %d\n", offset);
	  CDS_FSEEK( srcStore, (off_t) offset, SEEK_SET);  // offset should be 17038664 for frag 304262
	  
	  safeWrite( srcStore, fsread, 56);  
	  fclose( srcStore );
	}
	getFragStore( ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);

	fprintf( stderr, "set frag %d clr_bgn to %d, clr_end to %d\n",
			 fragIid, new_clr_bgn, new_clr_end);
  }

  return ( setStatus ); 
}

// this transfers the updated values from the nodes in the CIGraph to the containing contig
void updateIntUnitigPoss( NodeCGW_T *contig )
{
  int i, numUnitigs;
  MultiAlignT *cma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE);
  MultiAlignT *new_cma = CreateEmptyMultiAlignT();
  
  numUnitigs = GetNumIntUnitigPoss( cma->u_list );

  if ( UnitigPositions == NULL )
	UnitigPositions = CreateVA_IntElementPos(numUnitigs);
  ResetVA_IntElementPos( UnitigPositions );

  for ( i = 0; i < numUnitigs; i++)
  {
  	IntUnitigPos *upos = GetIntUnitigPos( cma->u_list, i);
	NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
	IntElementPos unitigPos;

	fprintf( stderr, "changing unitig %8d pos from (%10d, %10d) to (%10d, %10d)\n",
			 unitig->id, upos->position.bgn, upos->position.end, 
			 (int) unitig->offsetAEnd.mean, (int) unitig->offsetBEnd.mean);
	upos->position.bgn = unitig->offsetAEnd.mean;
	upos->position.end = unitig->offsetBEnd.mean;

	// load up UnitigPositions
	unitigPos.ident = unitig->id;
	unitigPos.type = AS_UNITIG;
	unitigPos.position.bgn = unitig->offsetAEnd.mean;
	unitigPos.position.end = unitig->offsetBEnd.mean;
	AppendIntElementPos( UnitigPositions, &unitigPos);
  }

  // ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, new_cma, contig->id, FALSE);

  new_cma = MergeMultiAlignsFast_new( ScaffoldGraph->sequenceDB,
									  ScaffoldGraph->fragStore,
									  UnitigPositions, FALSE, TRUE, GlobalData->aligner);

  // RemoveMultiAlignFromStore( ScaffoldGraph->sequenceDB, contig->id);
  // DeleteMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE);
  UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE);
  InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, contig->id, FALSE, new_cma, TRUE);
}

# if 0
void extendUnitigs( NodeCGW_T *unitig, int fragIid, extendableFrag extFrag, int leftContig)
{
  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd
  int unitigLeftEnd, unitigRightEnd, unitigOrientation;
  int contigLeftEnd, contigRightEnd, contigOrientation;
  
  // GetContigPositionInScaffold( contig, &contigLeftEnd, &contigRightEnd, &contigOrientation);
  GetCIPositionInScaffold( unitig, &unitigLeftEnd, &unitigRightEnd, &unitigOrientation);

  unitig->bpLength.mean += extFrag.extension;
  
  if ( leftContig == TRUE )
  {
	if ( unitigOrientation == 0 /* A_B */)
	{
	  fprintf (stderr, "extending AEnd of unitig %8d from %.0f to %.0f\n",
			   unitig->id, unitig->offsetAEnd.mean, unitig->offsetAEnd.mean + extFrag.extension);
	  unitig->offsetBEnd.mean += extFrag.extension;
	}
	else // unitigOrientation == 1 /* B_A */
	{
	  fprintf (stderr, "extending AEnd of unitig %8d from %.0f to %.0f\n",
			   unitig->id, unitig->offsetAEnd.mean, unitig->offsetAEnd.mean + extFrag.extension);
	  unitig->offsetAEnd.mean += extFrag.extension;
	}
  }
  // we're handling the right contig of the gap
  // have to step through all the unitigs of the right contig and add extFrag.extension to both offsets
  // except for the unitig of interest - one of it's ends stays at 0.0, the other gets added
  else 
  {
	ContigTIterator contigIter;
	NodeCGW_T *currUnitig;
	ContigT *contig = GetGraphNode( ScaffoldGraph->ContigGraph, unitig->info.CI.contigID);
	
	InitContigTIterator( ScaffoldGraph, contig->id, TRUE, FALSE, &contigIter);
	while ( NextContigTIterator( &contigIter ))
	{
	  currUnitig = GetGraphNode( ScaffoldGraph->CIGraph, contigIter.curr);

	  if (currUnitig->id != unitig->id)  // adjust everybody upwards except extended unitig
	  {
		currUnitig->offsetAEnd.mean += extFrag.extension;
		currUnitig->offsetBEnd.mean += extFrag.extension;
	  }
	  else
	  {
		if ( unitigOrientation == 0 /* A_B */)
		  currUnitig->offsetBEnd.mean += extFrag.extension;
		else  // unitig is B_A in contig
		  unitig->offsetAEnd.mean += extFrag.extension;
	  }
	}
  }
}
#endif

void DumpContigMultiAlignInfo ( int contigID )
{
#ifndef DIAG_PRINTS
  return;
#else
  
  fprintf( stderr, "------------------------------------------------------------\n");
  fprintf( stderr, "in DumpContigMultiAlignInfo, dumping info on contig %8d\n", contigID );
  
  cma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contigID, FALSE);

  fprintf( stderr, "in DumpContigMultiAlignInfo, contig %8d, strlen( consensus ): %9ld\n",
		   contigID, strlen( Getchar( cma->consensus, 0 )));

  for ( i = 0; i < GetNumIntMultiPoss( cma->f_list ); i++) {
       IntMultiPos *pos = GetIntMultiPos(cma->f_list,i);
       fprintf( stderr, "   in DumpContigMultiAlignInfo, fragment %8d, bgn: %10d, end: %10d, length: %10d\n", 
			pos->ident, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
  }
  for ( i = 0; i < GetNumIntUnitigPoss( cma->u_list ); i++)
  {
	IntUnitigPos *pos = GetIntUnitigPos( cma->u_list, i);
	NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, pos->ident);

	fprintf( stderr, "in DumpContigMultiAlignInfo, unitig %8d, bgn: %10d, end: %10d, length: %10d\n", 
			 unitig->id, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
	
	DumpUnitigMultiAlignInfo( unitig->id );	
  }
  fprintf( stderr, "\n");
  #endif
}

void DumpContigMultiAlignInfoDirect ( MultiAlignT *cma, int contigID )
{
  // MultiAlignT *cma = CreateEmptyMultiAlignT();

#ifndef DIAG_PRINTS
  return;
#else
  
  fprintf( stderr, "------------------------------------------------------------\n");
  fprintf( stderr, "in DumpContigMultiAlignInfoDirect, dumping info on contig %8d\n", contigID );
  
  // cma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contigID, FALSE);

  fprintf( stderr, "in DumpContigMultiAlignInfo, contig %8d, strlen( consensus ): %9ld\n",
		   contigID, strlen( Getchar( cma->consensus, 0 )));

  for ( i = 0; i < GetNumIntUnitigPoss( cma->u_list ); i++)
  {
	IntUnitigPos *pos = GetIntUnitigPos( cma->u_list, i);
	NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, pos->ident);

	fprintf( stderr, "in DumpContigMultiAlignInfo, unitig %8d, bgn: %10d, end: %10d, length: %10d\n", 
			 unitig->id, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
	
	DumpUnitigMultiAlignInfo( unitig->id );	
  }
  fprintf( stderr, "\n");
  #endif
}

void DumpUnitigMultiAlignInfo ( int unitigID )
{
#ifndef DIAG_PRINTS
  return;
#else
  
  fprintf( stderr, "in DumpUnitigMultiAlignInfo, dumping info on unitig %8d\n", unitigID );
  
  uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitigID, TRUE);

  fprintf( stderr, "in DumpUnitigMultiAlignInfo, unitig %8d, strlen( consensus ): %9ld\n",
		   unitigID, strlen( Getchar( uma->consensus, 0 )));

  for ( i = 0; i < GetNumIntMultiPoss( uma->f_list ); i++)
  {
	IntMultiPos *pos = GetIntMultiPos( uma->f_list, i);

	fprintf( stderr, "in DUMAI, fragment %8d, bgn: %10d, end: %10d, length: %10d, source: %d\n", 
			 pos->ident, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end),
			 pos->source);
  }
  #endif
}

static VA_TYPE(int32) *UngappedOffsets = NULL;

void DumpContigUngappedOffsets( int contigID )
{
  // int32 contigID = contig->id;


#ifndef DIAG_PRINTS
  return;
#else

  numCIs = GetNumIntUnitigPoss( ma->u_list );

  if( !UngappedOffsets )
  {
	UngappedOffsets = CreateVA_int32( 1000 );
  }
  //  fprintf(GlobalData->stderrc,"* UpdateNodeUnitigs for contig %d\n", contig->id);
  
  GetMultiAlignUngappedOffsets( ma, UngappedOffsets);
  offsets = Getint32( UngappedOffsets, 0);

  for( i = 0; i < numCIs ; i++)
  {
	IntUnitigPos *pos = GetIntUnitigPos( ma->u_list, i);
	NodeCGW_T *node = GetGraphNode( ScaffoldGraph->CIGraph, pos->ident);
	int flip = (pos->position.end < pos->position.bgn);
	int bgn, end;

	// mp->position is an interval.  We need to subtract one from
	// the upper end of the interval
	if( flip )
	{
	  bgn = pos->position.bgn - 1;
	  end = pos->position.end;
	}
	else
	{
	  bgn = pos->position.bgn;
	  end = pos->position.end - 1;
	}

	fprintf( stderr, "in DCUO, unitig %d, (bgn, end): ( %6d, %6d), offsets[bgn]: %10d, offsets[bgn]: %10d\n",
			 node->id, bgn, end, offsets[bgn], offsets[end] );
	fprintf( stderr, "in DCUO, unitig %d, pos->position.bgn: %10d, pos->position.end: %10d\n",
			 node->id, pos->position.bgn, pos->position.end );
  }

  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
  {
	IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
	int fragID = (int32)mp->source; // GetInfoByIID(ScaffoldGraph->iidToFragIndex, mp->ident)->fragIndex;
	// hack next line replaces above
	// int fragID = GetInfoByIID(ScaffoldGraph->iidToFragIndex, mp->ident)->fragIndex;
	CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragID);
	LengthT offset3p, offset5p;
	int ubgn, uend;
	int flip = (mp->position.end < mp->position.bgn);
	int bgn, end;
	
	// mp->position is an interval.  We need to subtract one from
	// the upper end of the interval
	if(flip)
	{
	  bgn = mp->position.bgn - 1;
	  end = mp->position.end;
	}
	else
	{
	  bgn = mp->position.bgn;
	  end = mp->position.end - 1;
	}

	fprintf( stderr, "in DCUO, contig %8d, frag %10d, mp->position.bgn: %10d, mp->position.end: %10d, "
			 "len: %10d, contained: %8d, source: %d\n",
			 contigID, frag->iid, mp->position.bgn, mp->position.end, abs(mp->position.bgn - mp->position.end),
			 mp->contained, (int32) mp->source);
  }
#endif
}

void getAlteredFragPositions( NodeCGW_T *unitig, fragPositions **fragPoss, int alteredFragIid, int extension)
{
  fragPositions *localFragPoss;
  int i, alteredFragIndex = NULLINDEX, orientation = 0, delta;
  MultiAlignT *uma = CreateEmptyMultiAlignT();

  // currently code does not handle negative extensions, ie, trimming
  if (extension <= 0)
	fprintf( stderr, "negative extension: %d\n", extension);
  
  // assert ( extension > 0 );
  
  uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE);

  localFragPoss = (fragPositions *) safe_malloc(  GetNumIntMultiPoss( uma->f_list ) * sizeof( fragPositions ));
  if ( localFragPoss == NULL)
  {
	fprintf( stderr, "failed to safe_malloc space for localFragPoss!\n");
	assert(0);
  }

  // get the current positions
  for ( i = 0; i < GetNumIntMultiPoss( uma->f_list ); i++)
  {
	IntMultiPos *pos = GetIntMultiPos( uma->f_list, i);
  
	localFragPoss[i].bgn = pos->position.bgn;
	localFragPoss[i].end = pos->position.end;

	if (pos->ident == alteredFragIid)
	{
	  alteredFragIndex = i;
 
	  if (pos->position.bgn < pos->position.end)
		orientation = 0;
	  else
		orientation = 1;
	}	
  }

  // alteredFragDelta is the amount by which the altered frag gets extended
  if (orientation == 0)  // nothing much to do, just extend the altered fragment
  {
	// first find out how much must be added to the frag to get to the end of the untig
	// int alteredFragDelta = unitig->bpLength.mean - localFragPoss[alteredFragIndex].end;
	int alteredFragDelta = strlen( Getchar( uma->consensus, 0)) - localFragPoss[alteredFragIndex].end;
	// then add how far it extends out into the gap
	alteredFragDelta += extension;
	localFragPoss[ alteredFragIndex ].end += alteredFragDelta;
  }
  else  // all frag positions get bumped by extension, that's how far the altered frag extends into the gap
  {
	// first find out how much must be added to the frag to get to the end of the untig
	// alteredFragDelta = localFragPoss[alteredFragIndex].end;
	// then add how far it extends into the gap
	// alteredFragDelta += extension;
	// extend the fragment
	// localFragPoss[ alteredFragIndex ].end -= alteredFragDelta;

	// this is the new minimum position in the unitig
	localFragPoss[ alteredFragIndex ].end = -extension;
	
	// if he extends off the front of the unitig, adjust everybody upward
	delta = localFragPoss[ alteredFragIndex ].end;
	if ( delta < 0 )
	{
	  for ( i = 0; i < GetNumIntMultiPoss( uma->f_list ); i++)
	  {
		localFragPoss[i].bgn -= delta;
		localFragPoss[i].end -= delta;
	  }
	}
  }
  *fragPoss = localFragPoss;
}

void getAlteredFragPositions_new( NodeCGW_T *unitig, fragPositions *fragPoss, int alteredFragIid, int extension)
{
  // fragPositions* localFragPoss = *fragPoss;
  int i, alteredFragIndex = NULLINDEX, orientation = 0, delta;
  MultiAlignT *uma = CreateEmptyMultiAlignT();
  // static int numFragments = -1;
  
  // currently code does not handle negative extensions, ie, trimming
  if (extension <= 0)
	fprintf( stderr, "negative extension: %d\n", extension);
  
  // assert ( extension > 0 );
  
  uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE);

  if ( GetNumIntMultiPoss( uma->f_list ) > FRAG_POSS_SIZE)
  {
	fprintf( stderr, "number of frags (%d) exceeds FRAG_POSS_SIZE (%d)\n", 
                 (int) GetNumIntMultiPoss( uma->f_list ), FRAG_POSS_SIZE);
	assert(0);
  }
  
  // get the current positions
  for ( i = 0; i < GetNumIntMultiPoss( uma->f_list ); i++)
  {
	IntMultiPos *pos = GetIntMultiPos( uma->f_list, i);
  
	fragPoss[i].bgn = pos->position.bgn;
	fragPoss[i].end = pos->position.end;

	if (pos->ident == alteredFragIid)
	{
	  alteredFragIndex = i;
 
	  if (pos->position.bgn < pos->position.end)
		orientation = 0;
	  else
		orientation = 1;
	}	
  }

  // alteredFragDelta is the amount by which the altered frag gets extended
  if (orientation == 0)  // nothing much to do, just extend the altered fragment
  {
	// first find out how much must be added to the frag to get to the end of the untig
	// int alteredFragDelta = unitig->bpLength.mean - fragPoss[alteredFragIndex].end;
	int alteredFragDelta = strlen( Getchar( uma->consensus, 0)) - fragPoss[alteredFragIndex].end;
	// then add how far it extends out into the gap
	alteredFragDelta += extension;
	fragPoss[ alteredFragIndex ].end += alteredFragDelta;
  }
  else  // all frag positions get bumped by extension, that's how far the altered frag extends into the gap
  {
	// first find out how much must be added to the frag to get to the end of the untig
	// alteredFragDelta = localFragPoss[alteredFragIndex].end;
	// then add how far it extends into the gap
	// alteredFragDelta += extension;
	// extend the fragment
	// localFragPoss[ alteredFragIndex ].end -= alteredFragDelta;

	// this is the new minimum position in the unitig
	fragPoss[ alteredFragIndex ].end = -extension;
	
	// if he extends off the front of the unitig, adjust everybody upward
	delta = fragPoss[ alteredFragIndex ].end;
	if ( delta < 0 )
	{
	  for ( i = 0; i < GetNumIntMultiPoss( uma->f_list ); i++)
	  {
		fragPoss[i].bgn -= delta;
		fragPoss[i].end -= delta;
	  }
	}
  }
}

// this is fast, things are mostly in order
void bubbleSortIUMs( IntMultiPos *f_list, int numIMPs)
{
  IntMultiPos tempIMP;
  int i, j;
  
  for ( i = numIMPs; --i > 0; )
  {
	for ( j = 0; j < i; j++)
	{	  
	  if ( min( f_list[j].position.bgn, f_list[j].position.end ) > 
		   min( f_list[j+1].position.bgn, f_list[j+1].position.end ))
	  {
		if ( f_list[0].ident == 300268 )
		  fprintf( stderr, "swapping %2d (frag %8d) and %2d (frag %8d)\n", i, f_list[i].ident, j, f_list[j].ident);

		memcpy( &tempIMP, &f_list[j], sizeof( IntMultiPos ) );
		memcpy( &f_list[j], &f_list[j+1], sizeof( IntMultiPos ) );
		memcpy( &f_list[j+1], &tempIMP, sizeof( IntMultiPos ) );
	  }
	}
  }
}

void SynchUnitigTWithMultiAlignT( NodeCGW_T *unitig )
{
  // IntUnitigPos *pos = GetIntUnitigPos( ma->u_list, i);
  // NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, pos->ident);
  MultiAlignT *uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE);
  
  fprintf( stderr, "in SynchUnitigTWithMultiAlignT, for unitig %d strlen( ma->consensus ) = " F_SIZE_T "\n",
		   unitig->id, strlen( Getchar( uma->consensus, 0) ));
  
  // unitig->bpLength.mean = strlen( Getchar( uma->consensus, 0) );      // set length
  unitig->bpLength.mean = GetMultiAlignUngappedLength(uma);
#if 0
  if ( unitig->offsetAEnd.mean < unitig->offsetBEnd.mean )          // ordering info is okay
  {
	unitig->offsetAEnd.mean = pos->position.bgn;
	unitig->offsetBEnd.mean = pos->position.end;
  }
  else 
  {
	unitig->offsetAEnd.mean = pos->position.end;
	unitig->offsetBEnd.mean = pos->position.bgn;
  }
#endif  
}

// this routine shifts the frag with iid extendedFragIid to the front of the f_list
void leftShiftIUM( IntMultiPos *f_list, int numFrags, int extendedFragIid)
{
  int i, currPos = 0, numShiftedInUnitig = 0;
  IntMultiPos tempIMP;
  
  // find out the extended frag's current position
  for ( i = 0; i < numFrags; i++)
	if ( f_list[i].ident == extendedFragIid )
	  currPos = i;
  
  // save the extended frag off to the side
  memcpy( &tempIMP, &f_list[currPos], sizeof( IntMultiPos ) );

  // now shift everybody below currPos up a position
  for ( i = currPos - 1; i >= 0; i--)
  {
	memcpy( &f_list[i+1], &f_list[i], sizeof( IntMultiPos ) );
	numShiftedInUnitig++;
  }
  
  // and place tempPos in position 0 of the unitig
  memcpy( &f_list[0], &tempIMP, sizeof( IntMultiPos ) );
  
#if 0
  // now do the same shifting for the contig
  {
	NodeCGW_T *contig = GetGraphNode( ScaffoldGraph->ContigGraph, contigID);
	MultiAlignT *cma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contigID, FALSE);
	IntMultiPos *cf_list = cma->f_list;
	int icnt, numContigFrags = GetNumIntMultiPoss( cma->f_list );
	
	// find out the extended frag's current position
	for ( i = 0; i < numContigFrags; i++)
	  if ( cf_list[i].ident == extendedFragIid )
		currPos = i;
	
	// save the extended frag off to the side
	memcpy( &tempIMP, &cf_list[currPos], sizeof( IntMultiPos ) );
	
	// now shift the appropriate number of frags up a position starting at currPos
	for ( i = currPos - 1, icnt = 0; icnt < numShiftedInUnitig; i--, icnt++)
	{
	  memcpy( &cf_list[i+1], &cf_list[i], sizeof( IntMultiPos ) );
	  numShiftedInUnitig++;
	}
	
	// and place tempPos in position 0
	memcpy( &cf_list[currPos - numShiftedInUnitig], &tempIMP, sizeof( IntMultiPos ) );
  }
#endif
}

// this routine shifts the frag with iid extendedFragIid to the appropriate place in the f_list
void rightShiftIUM( IntMultiPos *f_list, int numFrags, int extendedFragIid)
{
  int i, currPos = 0;
  IntMultiPos tempIMP;
  int numShifted, shiftedFrag, fragToMovePos = NULLINDEX, numToShift;
  
  // find out the extended frag's current position
  for ( i = 0; i < numFrags; i++)
	if ( f_list[i].ident == extendedFragIid )
	  currPos = i;  

  // we need to move all frags that 
  // i) start to the left of the extended frag and 
  // ii) are not contained by another frag 
  // before the extended frag in the array.  We could move the contained ones, but don't have to

  numShifted = 0;
  shiftedFrag = TRUE;
  while ( shiftedFrag )
  {
	shiftedFrag = FALSE;
	
	// look for a candidate frag
	i = currPos + numShifted + 1;
	while ( i < numFrags )
	{
	  if ( min( f_list[i].position.bgn, f_list[i].position.end ) < 
		   min(f_list[currPos + numShifted].position.bgn, f_list[currPos + numShifted].position.end)
			&& f_list[i].contained == FALSE)
	  {
		fragToMovePos = i;
		shiftedFrag = TRUE;
		break;
	  }
	  i++;
	}

	if (shiftedFrag)
	{
	  // save the extended frag to be shifted off to the side
	  memcpy( &tempIMP, &f_list[fragToMovePos], sizeof( IntMultiPos ) );
	  
	  // now shift all the frags between currPos + numShifted up a position in the array
	  numToShift = fragToMovePos - ( currPos + numShifted );
	  i = 0;
	  while ( i < numToShift)
	  {
		memcpy( &f_list[fragToMovePos - i], &f_list[fragToMovePos - i - 1], sizeof( IntMultiPos ) );
		i++;
	  }
	  // and tempPos into open position
	  memcpy( &f_list[currPos + numShifted], &tempIMP, sizeof( IntMultiPos ) );

	  numShifted++;
	}
  }
}

static MultiAlignT *savedLeftUnitigMA = NULL;
static MultiAlignT *savedRightUnitigMA = NULL;
static MultiAlignT *savedLeftContigMA = NULL;
static MultiAlignT *savedRightContigMA = NULL;

void saveFragAndUnitigData( int lFragIid, int rFragIid )
{
  CIFragT *leftFrag, *rightFrag;
  InfoByIID *info;
  
  if (savedLeftUnitigMA == NULL)
  {
	savedLeftUnitigMA = CreateEmptyMultiAlignT();
	savedRightUnitigMA = CreateEmptyMultiAlignT();
	savedLeftContigMA = CreateEmptyMultiAlignT();
	savedRightContigMA = CreateEmptyMultiAlignT();
  }

  // grab the multialignments
  if ( lFragIid != -1)
  {
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
	assert(info->set);
	leftFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);

	ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, savedLeftUnitigMA, leftFrag->cid, TRUE);
	ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, savedLeftContigMA, leftFrag->contigID, FALSE);
  }
  
  if ( rFragIid != -1)
  {
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
	assert(info->set);
	rightFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);

	ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, savedRightUnitigMA, rightFrag->cid, TRUE);
	ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, savedRightContigMA, rightFrag->contigID, FALSE);
  }
}

void restoreFragAndUnitigData( int lFragIid, int rFragIid )
{
  NodeCGW_T *unitig;
  CIFragT *leftFrag, *rightFrag;
  
  // first the left frag
  if ( lFragIid != -1)
  {
	InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);

	assert(info->set);
	leftFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
	unitig = GetGraphNode( ScaffoldGraph->CIGraph, leftFrag->cid);
	UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE);
	InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE, savedLeftUnitigMA, FALSE);
	SynchUnitigTWithMultiAlignT( unitig );

	UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, leftFrag->contigID, FALSE);
	InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, leftFrag->contigID, FALSE, savedLeftContigMA, FALSE);

	fprintf( stderr, "in restoreFragAndUnitigData, left contig %d:\n", leftFrag->contigID);
	DumpContigMultiAlignInfo ( leftFrag->contigID );
	DumpContigUngappedOffsets( leftFrag->contigID );
  }
  
  // now the right frag
  if ( rFragIid != -1)
  {
	InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);

	assert(info->set);
	rightFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
	unitig = GetGraphNode( ScaffoldGraph->CIGraph, rightFrag->cid);
	UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE);
	InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, unitig->id, TRUE, savedRightUnitigMA, FALSE);
	SynchUnitigTWithMultiAlignT( unitig );
	
	UnloadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, rightFrag->contigID, FALSE);
	InsertMultiAlignTInSequenceDB( ScaffoldGraph->sequenceDB, rightFrag->contigID, FALSE, savedRightContigMA, FALSE);

	fprintf( stderr, "in restoreFragAndUnitigData, right contig %d:\n", rightFrag->contigID);
	DumpContigMultiAlignInfo ( rightFrag->contigID );
	DumpContigUngappedOffsets( rightFrag->contigID );
  }
}

static int MaxBegGapSaved, MaxEndGapSaved, MaxGapsSaved, MaxInteriorGapSaved, asymmetricEndsSaved;
static int DefaultMaxBegGap, DefaultMaxEndGap, DefaultMaxGaps, DefaultMaxInteriorGap, DefaultAsymmetricEnds;

void saveDefaultLocalAlignerVariables( void )
{
  DefaultMaxBegGap = MaxBegGap;
  DefaultMaxEndGap = MaxEndGap;
  DefaultMaxGaps = MaxGaps;
  DefaultMaxInteriorGap = MaxInteriorGap;
  DefaultAsymmetricEnds = asymmetricEnds;
}

void saveLocalAlignerVariables( void )
{
  MaxBegGapSaved = MaxBegGap;
  MaxEndGapSaved = MaxEndGap;
  MaxGapsSaved = MaxGaps;
  MaxInteriorGapSaved = MaxInteriorGap;
  asymmetricEndsSaved = asymmetricEnds;
  
  MaxBegGap = DefaultMaxBegGap;
  MaxEndGap = DefaultMaxEndGap;
  MaxGaps = DefaultMaxGaps;
  MaxInteriorGap = DefaultMaxInteriorGap;
  asymmetricEnds = DefaultAsymmetricEnds;
}

void restoreLocalAlignerVariables( void )
{
  MaxBegGap = MaxBegGapSaved;
  MaxEndGap = MaxEndGapSaved;
  MaxGaps = MaxGapsSaved;
  MaxInteriorGap = MaxInteriorGapSaved;
}



#if 0
					if (0)
					{
					  // adjust the iterator to reflect the new contig
					  InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
					  
					  // go until we hit the node we just inserted
					  fprintf( stderr, "during gap closing of scaffold %d\n", scaff->id);
					  while ( NextCIScaffoldTIterator( &CIsTemp ) && 
							  CIsTemp.next != GetNumGraphNodes( ScaffoldGraph->ContigGraph ) - 1)
					  {
						fprintf( stderr, "prev, curr, next: %d, %d, %d\n", CIsTemp.prev, CIsTemp.curr, CIsTemp.next);
					  }
					  fprintf( stderr, "NumContigs: %d\n", GetNumGraphNodes( ScaffoldGraph->ContigGraph ));
					  fprintf( stderr, "CIsTemp.prev: %d\n", CIsTemp.prev);
					  fprintf( stderr, "CIsTemp.curr: %d\n", CIsTemp.curr);
					  fprintf( stderr, "CIsTemp.next: %d\n", CIsTemp.next);					  
					}
					
#endif

#if 0
void resetSourceFields( MultiAlignT *ma )
{
  int icnt, numNodeFrags = GetNumIntMultiPoss( ma->f_list );

  for ( icnt = 0; icnt < numNodeFrags; icnt++)
  {
	IntMultiPos *pos = GetIntMultiPos( ma->f_list, icnt);
	int fragIndex = GetInfoByIID( ScaffoldGraph->iidToFragIndex, pos->ident)->fragIndex;

	pos->source = (char *)fragIndex; // CMM: this looks weird.
  }
}
#endif

void printGapSizes()
{
  int sid;
  
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
	CIScaffoldTIterator CIsTemp;
	CIScaffoldT * scaff;
	int icnt;
	LengthT gapEstimate;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    
	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds    
	if ((isDeadCIScaffoldT(scaff)) ||
		(scaff->type != REAL_SCAFFOLD) ||
		(scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}
	fprintf(stderr,"\n=====================================================================\n");
	fprintf(stderr,"=== examining scaffold %d, size %f\n", sid, scaff->bpLength.mean);
	  
	icnt = 0;
    InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,
							FALSE, &CIsTemp);
    while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ChunkInstanceT *lchunk, *rchunk;

	  // not walking off of scaffolds currently
	  if (CIsTemp.next == -1)
		break;
	  
      // find the chunks in the gap by walking between the chunk <CIs.curr>
      // and the chunk <CIs.next>
      lchunk = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
      assert(lchunk != NULL);

	  rchunk = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.next);
      assert(rchunk != NULL);

	  gapEstimate = FindGapLength(lchunk, rchunk, TRUE);

	  fprintf( stderr, "gap between contigs %8d and %8d is %8.0f\n",
			   lchunk->id, rchunk->id, gapEstimate.mean);
	}
  }
}

int writeEcrCheckpoint( int *numGapsInScaffold, int *numGapsClosedInScaffold,
						int *numSmallGapsInScaffold, int *numSmallGapsClosedInScaffold,
						int *numLargeGapsInScaffold, int *numLargeGapsClosedInScaffold)
{
  char ckpFileName[1024];
  FILE *ckpFile;
  int i;
  
  sprintf( ckpFileName, "%s.ecr.ckp.%d", GlobalData->File_Name_Prefix, ScaffoldGraph->checkPointIteration - 1);

  ckpFile = fopen( ckpFileName, "w+");
  if (ckpFile == NULL)
  {
	fprintf( stderr, "Could not open checkpoint file %s, aborting.\n", ckpFileName);
	assert(0);
  }
  
  fprintf( ckpFile, "%d\n", GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph ));
  for ( i = 0; i < GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph ); i++)
  {
	CIScaffoldT *scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, i);
	assert( scaff != NULL );
	
	//if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD) ||	(scaff->info.Scaffold.numElements < 2))
	//continue;
	
	fprintf( ckpFile, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", scaff->id, 
			 numGapsInScaffold[ scaff->id ], numGapsClosedInScaffold[ scaff->id ],
			 numSmallGapsInScaffold[ scaff->id ], numSmallGapsClosedInScaffold[ scaff->id ],
			 numLargeGapsInScaffold[ scaff->id ], numLargeGapsClosedInScaffold[ scaff->id ] );
  }
  fclose( ckpFile );
  return 0;
}

int loadEcrCheckpoint( int ckptNum, int *numGapsInScaffold, int *numGapsClosedInScaffold,
					   int *numSmallGapsInScaffold, int *numSmallGapsClosedInScaffold,
					   int *numLargeGapsInScaffold, int *numLargeGapsClosedInScaffold)
{
  char ckpFileName[1024];
  FILE *ckpFile;
  int i, scaffid;
  int totalGaps = 0, totalGapsClosed = 0, 
	totalSmallGaps = 0, totalSmallGapsClosed = 0, totalLargeGaps = 0, totalLargeGapsClosed = 0;
  int numGraphNodes;
  
  sprintf( ckpFileName, "%s.ecr.ckp.%d", GlobalData->File_Name_Prefix, ckptNum);

  ckpFile = fopen( ckpFileName, "r");
  if (ckpFile == NULL)
  {
	fprintf( stderr, "Warning: could not open checkpoint file %s for reading, continuing.\n", ckpFileName);
	return -1;
  }
  
  fscanf( ckpFile, "%d\n", &numGraphNodes);
  for ( i = 0; i < GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph ); i++)
  {
	// CIScaffoldT* scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, i);
	// assert( scaff );
	//if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD) ||	(scaff->info.Scaffold.numElements < 2))
	//continue;
	
	fscanf( ckpFile, "%d %d %d %d %d %d %d\n", &scaffid, 
			 &numGapsInScaffold[ i ], &numGapsClosedInScaffold[ i ],
			 &numSmallGapsInScaffold[ i ], &numSmallGapsClosedInScaffold[ i ],
			 &numLargeGapsInScaffold[ i ], &numLargeGapsClosedInScaffold[ i ] );

	//if ( scaffid != i)
	 // fprintf( stderr, "Warning: scaffid (%d) != i (%d) in loadEcrCheckpoint!\n",
	//		   scaffid, i);

	totalGaps += numGapsInScaffold[ scaffid ];
	totalGapsClosed += numGapsClosedInScaffold[ scaffid ];
	totalSmallGaps += numSmallGapsInScaffold[ scaffid ];
	totalSmallGapsClosed += numSmallGapsClosedInScaffold[ scaffid ];
	totalLargeGaps += numLargeGapsInScaffold[ scaffid ];
	totalLargeGapsClosed += numLargeGapsClosedInScaffold[ scaffid ];
  }
  fclose( ckpFile );

  fprintf( stderr, "Stats after loadEcrCheckpoint for ckptNum %d:\n", ckptNum);
  fprintf( stderr, "           totalGaps: %d\n", totalGaps);
  fprintf( stderr, "     totalGapsClosed: %d\n", totalGapsClosed);
  fprintf( stderr, "      totalSmallGaps: %d\n", totalSmallGaps);
  fprintf( stderr, "totalSmallGapsClosed: %d\n", totalSmallGapsClosed);
  fprintf( stderr, "      totalLargeGaps: %d\n", totalLargeGaps);
  fprintf( stderr, "totalLargeGapsClosed: %d\n", totalLargeGapsClosed);
  
  return 0;
}

