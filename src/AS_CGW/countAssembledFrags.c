
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
static char CM_ID[] = "$Id: countAssembledFrags.c,v 1.1.1.1 2004-04-14 13:50:32 catmandew Exp $";


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
#include "Utils_REZ.h"
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
#include "Stats_CGW.h"
#include "GapWalkerREZ.h"
#include "FbacREZ.h"
#include "countAssembledFrags.h"

int compFrags( const void *s1, const void *s2)
{
  const int * t1 = s1;
  const int * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( *t1 < *t2)
    return -1;
  else if ( *t1 > *t2)
    return 1;
  else 
    return 0;
}

int compRanges( const void *s1, const void *s2)
{
  const scaffoldRangeT * t1 = s1;
  const scaffoldRangeT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( t1->firstRange < t2->firstRange )
    return -1;
  else if ( t1->firstRange > t2->firstRange )
    return 1;
  else 
  {
    assert(0);   // ranges are unique
    return 0;
  }
}

int compContigIDs( const void *s1, const void *s2)
{
  const fragDataT * t1 = s1;
  const fragDataT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( t1->contigID < t2->contigID )
    return -1;
  else   if ( t1->contigID > t2->contigID )
    return 1;
  else // in same contig, sort by contigOffset5p
  {
    if (t1->contigOffset5p < t2->contigOffset5p)
      return -1;
    else if (t1->contigOffset5p > t2->contigOffset5p)
      return 1;
    else
      return 0;
  }
}

int compScaffoldIids( const void *s1, const void *s2)
{
  const scaffoldIidRangeT * t1 = s1;
  const scaffoldIidRangeT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( t1->iidStart < t2->iidStart )
    return -1;
  else   if ( t1->iidStart > t2->iidStart )
    return 1;
  else 
    return 0;
}

int compConsistentRangesByScaffold( const void *s1, const void *s2)
{
  const rangeDataT * t1 = s1;
  const rangeDataT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if ( t1->scaffoldID < t2->scaffoldID )
    return -1;
  else if ( t1->scaffoldID > t2->scaffoldID )
    return 1;
  else 
  {
    if ( t1->scaffoldOffsetMin < t2->scaffoldOffsetMin )
      return -1;
    else if ( t1->scaffoldOffsetMin > t2->scaffoldOffsetMin )
      return 1;
    else
      return 0;
  }
}

int main(int argc, char *argv[])
{
  int32 restartFromCheckpoint = NULLINDEX;
  Global_CGW *data;
  char *inputPath;
  char *prefix;
  MesgReader reader;
  MesgWriter writer;
  MesgWriter errorWriter;
  FILE *myerr = stderr; 
  FILE *myout = stdout; 
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  int i, index;
  CIFragT *frag, *previousFrag;
  NodeCGW_T *unitig;
  ContigT *contig;
  CIScaffoldT *scaffold;
  int ifrag;
  long int numBACFrags = 0;
  long int numBACFragsInAss = 0;
  long int numBACFragsInWalkSurrogates = 0;
  long int numBACFragsInSurrogates = 0;
  long int numBACFragsNotInSurrogates = 0;
  long int numBACFragsInUniques = 0;
  double unitigLengthsBFNIS = 0.0;
  int numUnitigsBFNIS = 0;
  int *fragsNotInAss, prior, inRange, onEnd;
  CDS_CID_t gapSizeInFrags, beginFragIid, endFragIid;
  char *seenUnitig;
  int *localeFragsInUnitig, *localeFragsInContig, *localeFragsInScaffold;
  int *gapsFlanked;
  int flankingFragsInSameScaffold = 0;
  int numMissingBACFragRanges = 0,
    *numMissingBACFragRangesPerLocale;
  int flankingContigsAdjacent = 0;
  CDS_CID_t maxLocale;
  CDS_CID_t previousFragIid, previousUnitigId, previousContigId, previousScaffoldId, previousLocale;
  int numUnitigContainedGaps = 0, numContigContainedGaps = 0, numPureContigContainedGaps = 0;
  int *numFragsInLocale, *missingFragsPerLocale, *contigsPerLocale, *scaffoldGapsPerLocale, *fragGapsPerLocale;
  int *numUnitigContainedGapsPerLocale, // the number of missing frags whose flanking frags are in the same unitig
    *numContigContainedGapsPerLocale;   // the number of missing frags whose flanking frags are in the same contig
  CDS_COORD_t fragsScaffoldSeparation;
  int *sizesUnitigContainedGaps, *sizesContigContainedGaps;
  int numContigsInLocale = 0, numUnitigsInLocale = 0, numScaffoldsInLocale = 0, currLocale = -1;
  int missingFragsRecount = 0;
  int fragsMisoriented;
  
  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("optimizeDistLibs.stderr","w");
  
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:f:g:n:")) != EOF)){
      switch(ch) {
        case 'c':
        {
          strcpy( data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = 1;
          
        }
        break;
        case 'f':
        {
          strcpy( data->Frag_Store_Name, argv[optind - 1]);
          setFragStore = 1;
        }
        break;
        case 'g':
        {
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = 1;
        }
        break;	  
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }
    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0))
    {
      fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d outputPath = %s\n",
              argc, optind, setFragStore,setGatekeeperStore, outputPath);
      fprintf (stderr, "USAGE:  loadcgw -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n");
      exit (EXIT_FAILURE);
    }
  }
  
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);
  
  seenUnitig = (char *) safe_calloc ( GetNumGraphNodes(ScaffoldGraph->CIGraph), sizeof( char ));
  gapsFlanked = (int *) safe_calloc ( GetNumGraphNodes(ScaffoldGraph->ContigGraph), sizeof( int ));
  
  fprintf( stderr, "GetNumVA_CIFragT( ScaffoldGraph->CIFrags ): %d\n",
           (int) GetNumVA_CIFragT( ScaffoldGraph->CIFrags ));
  
  maxLocale = NULLINDEX;
  for (ifrag = 0; ifrag < GetNumVA_CIFragT( ScaffoldGraph->CIFrags ); ifrag++)
  {
    frag = GetCIFragT( ScaffoldGraph->CIFrags, ifrag);
    if (frag->locale >= 0)
    {
      numBACFrags++;
      unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
      contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
      
      if (contig->scaffoldID != NULLINDEX)
        numBACFragsInAss++;
      else if (unitig->info.CI.numInstances > 0)
      {
        numBACFragsInSurrogates += 1; // unitig->info.CI.numInstances;
      }
      else
      {
        numBACFragsNotInSurrogates += 1;
        
        if ( !seenUnitig[ unitig->id ])
        {
          unitigLengthsBFNIS += unitig->bpLength.mean;
          seenUnitig[ unitig->id ] = TRUE;
          numUnitigsBFNIS++;
        }
      }
      
      if (ifrag % 100 == 0 && 0)
        fprintf( stderr, "frag->contigID: " F_CID ", contig->scaffoldID: " F_CID "\n", 
                 frag->contigID, contig->scaffoldID);
      
      if (frag->locale >= maxLocale)
        maxLocale = frag->locale;
    }
  }
  
  fragsNotInAss = (int *) safe_malloc( numBACFragsNotInSurrogates * sizeof( int ));
  sizesUnitigContainedGaps = (int *) safe_malloc( numBACFragsNotInSurrogates * sizeof( int ));
  sizesContigContainedGaps = (int *) safe_malloc( numBACFragsNotInSurrogates * sizeof( int ));
  
  numFragsInLocale = (int *) safe_malloc( maxLocale * sizeof( int ));
  missingFragsPerLocale = (int *) safe_malloc( maxLocale * sizeof( int ));
  contigsPerLocale = (int *) safe_malloc( maxLocale * sizeof( int ));
  scaffoldGapsPerLocale = (int *) safe_malloc( maxLocale * sizeof( int ));
  numMissingBACFragRangesPerLocale = (int *) safe_malloc( maxLocale * sizeof( int ));
  numUnitigContainedGapsPerLocale = (int *) safe_malloc( maxLocale * sizeof( int ));  // unitig contains range of missing frags
  numContigContainedGapsPerLocale = (int *) safe_malloc( maxLocale * sizeof( int ));  // contig contains range of missing frags
  
  for (i = 1; i <= maxLocale; i++)
  {
    numFragsInLocale[i] = 0;
    missingFragsPerLocale[i] = 0;
    // these arrays should be filled when we walk across scaffolds, not here where we walk frags by index
    // contigsPerLocale[i] = 0;
    // scaffoldGapsPerLocale[i] = 0;
    // fragGapsPerLocale[i] = 0;
    numUnitigContainedGapsPerLocale[i] = 0;
    numContigContainedGapsPerLocale[i] = 0;	
    numMissingBACFragRangesPerLocale[i] = 0;
  }
  
  index = 0;
  for (ifrag = 0; ifrag < GetNumVA_CIFragT( ScaffoldGraph->CIFrags ); ifrag++)
  {
    frag = GetCIFragT( ScaffoldGraph->CIFrags, ifrag);
    if (frag->locale >= 0)
    {
      unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
      contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
      
      numFragsInLocale[ frag->locale ]++;
      
      // if frag is not in a scaffold or a surrogate he's missing from the assembly
      if ( !(contig->scaffoldID != NULLINDEX || unitig->info.CI.numInstances > 0))
      {
        fragsNotInAss[ index ] = frag->iid;
        missingFragsPerLocale[ frag->locale ]++;
        index++;
      }
    }
  }
  
  qsort( fragsNotInAss, numBACFragsNotInSurrogates, sizeof( int ), &compFrags);
  
  fprintf( stderr, "\nRanges of missing frags ************************************\n");
  
  prior = -2;
  for ( ifrag = 0; ifrag < numBACFragsNotInSurrogates; ifrag++)
  {
    InfoByIID *info;
    // fprintf( stderr, "frag: %d\n", fragsNotInAss[ ifrag ]);
    
    if ( fragsNotInAss[ ifrag ] != prior + 1)
    {
      int dec = -1;
      int inScaffold = FALSE;
      
      // get info on the frags that preceed the current frag
      // looking for the first one that's in a scaffold
      while (!inScaffold)
      {
        info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, fragsNotInAss[ ifrag ] + dec);
        assert(info->set);
        previousFrag = frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
        unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
        contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
        if (contig->scaffoldID != NULLINDEX)
        {
          inScaffold = TRUE;
          scaffold = GetGraphNode( ScaffoldGraph->ScaffoldGraph,
                                   contig->scaffoldID);
          if (scaffold->info.Scaffold.AEndCI == contig->id ||
              scaffold->info.Scaffold.BEndCI == contig->id)
            onEnd = TRUE;
          else
            onEnd = FALSE;
          fprintf( stderr, "          fragiid\t utg_id\t utg_len\t ctg_id\t scf_id\t onEnd\n" );
          fprintf( stderr, "previous: " F_CID "\t " F_CID "\t %lf\t " F_CID "\t " F_CID "\t %d\n",
                   frag->iid, unitig->id, unitig->bpLength.mean,
                   frag->contigID, contig->scaffoldID,
                   onEnd);
          previousFragIid = frag->iid;
          previousUnitigId = unitig->id;
          previousContigId = contig->id;
          previousScaffoldId = contig->scaffoldID;
          previousLocale = frag->locale;
        }
        else
          dec--;
      }
      
      // get info on the current frag
      info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, fragsNotInAss[ ifrag ]);
      assert(info->set);
      frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
      unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
      beginFragIid = frag->iid;
      numMissingBACFragRanges++;
      numMissingBACFragRangesPerLocale[ frag->locale ]++;
      fprintf( stderr, "   begin: " F_CID "\t " F_CID "\t %lf\t " F_CID "\t %d\n",
               frag->iid, unitig->id, unitig->bpLength.mean,
               frag->contigID, -1);
    }	  
    
    if ( (fragsNotInAss[ ifrag ] != fragsNotInAss[ ifrag + 1 ] - 1) || 
         (ifrag == numBACFragsNotInSurrogates - 1) )
    {
      // get info on this frag
      info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, fragsNotInAss[ ifrag ]);
      assert(info->set);
      frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
      unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
      endFragIid = frag->iid;
      fprintf( stderr, "  ending: " F_CID "\t " F_CID "\t %lf\t " F_CID "\t %d\n",
               frag->iid, unitig->id, unitig->bpLength.mean, 
               frag->contigID, -1);
      
      // get info on the frags that follows this one, providing it exists
      if (fragsNotInAss[ ifrag ] + 1 <= GetNumVA_CIFragT( ScaffoldGraph->CIFrags ))
      {
        int inc = 1;
        int inScaffold = FALSE;
        
        // get info on the frags that follow the current frag
        // looking for the first one that's in a scaffold
        while (!inScaffold && 
               fragsNotInAss[ ifrag ] + inc + 1 < GetNumVA_CIFragT( ScaffoldGraph->CIFrags ))
        {
          info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, fragsNotInAss[ ifrag ] + inc);
          assert(info->set);
          frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
          unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
          contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
          if (contig->scaffoldID != NULLINDEX)
          {
            inScaffold = TRUE;
            scaffold = GetGraphNode( ScaffoldGraph->ScaffoldGraph, contig->scaffoldID);
            if (scaffold->info.Scaffold.AEndCI == contig->id ||
                scaffold->info.Scaffold.BEndCI == contig->id)
              onEnd = TRUE;
            else
              onEnd = FALSE;
            fprintf( stderr, "    next: " F_CID "\t " F_CID "\t %lf\t " F_CID "\t " F_CID "\t %d\n",
                     frag->iid, unitig->id, unitig->bpLength.mean, 
                     frag->contigID, contig->scaffoldID,
                     onEnd);
            if (contig->scaffoldID == previousScaffoldId && 
                frag->contigID != previousContigId &&
                frag->locale == previousLocale)
            {
              fprintf( stderr, "flankingFragsInSameScaffold!\n");
              fprintf( stderr, "previousContigId: " F_CID ", contig->AEndNext: " F_CID ", contig->BEndNext: " F_CID ", frag->locale: " F_CID "\n",
                       previousContigId, contig->AEndNext, contig->BEndNext, frag->locale);
              if (contig->BEndNext == previousContigId)
                flankingContigsAdjacent++;
              
              // determine separation in scaffold
              // if (getScaffoldSeparation( previousFragIid, frag->iid, &fragsScaffoldSeparation, TRUE) != -1)
              {
                LengthT scaffoldGapSize;
                ContigT *previousContig = GetGraphNode( ScaffoldGraph->ContigGraph, previousContigId);
                
                // since we are moving along in fragIid space, not on scaffolds
                // we're not sure of relative postions of contig and previousContig in their scaffold
                if (contig->offsetAEnd.mean < previousContig->offsetAEnd.mean)
                  scaffoldGapSize = FindGapLength( contig, previousContig, FALSE);
                else
                  scaffoldGapSize = FindGapLength( previousContig, contig, FALSE);
                
                if (scaffoldGapSize.mean == -20)
                {
                  // fprintf( stderr, "scaffoldGapSize.mean: %f, scaffoldGapSize.variance: %f\n",
                  //   scaffoldGapSize.mean, scaffoldGapSize.variance);				
                  scaffoldGapSize.variance = 10.0;
                }
                
                // if (scaffoldGapSize.mean != -20)
                {
                  gapsFlanked[ flankingFragsInSameScaffold++ ] = (int) 
                    (abs( fragsScaffoldSeparation -
                          (frag->localePos.bgn - previousFrag->localePos.end)) / 
                     sqrt(scaffoldGapSize.variance));
                }
                
                fprintf( stderr, " fragsScaffoldSeparation: " F_COORD "\n", fragsScaffoldSeparation);
                fprintf( stderr, "   fragsLocaleSeparation: " F_COORD "\n", frag->localePos.bgn - previousFrag->localePos.end);
                fprintf( stderr, "    scaffoldGapSize.mean: %.2lf\n", scaffoldGapSize.mean);
                fprintf( stderr, "scaffoldGapSize.variance: %.2lf (%.2lf)\n", scaffoldGapSize.variance,
                         sqrt( scaffoldGapSize.variance ));
              }
            }
          }
          else
            inc++;
        }
      }
      
      gapSizeInFrags = endFragIid - beginFragIid + 1;
      fprintf( stderr, "gap size: " F_CID " (" F_CID "," F_CID ")\n",
               gapSizeInFrags, beginFragIid, endFragIid);
      // if (gapSizeInFrags > 1)
      //	largegap = TRUE;
      
      fragsMisoriented = FALSE;
      if ( previousContigId == contig->id && previousLocale == frag->locale)
      {
        fprintf( stderr, "***--- contig contained gap ---*** (" F_CID "," F_CID ")\n", beginFragIid, endFragIid);
        if (previousLocale == frag->locale)
        {
          int fragsScaffoldSeparation;
          
          // a -1 return from getScaffoldSeparation indicates problems
          if (getScaffoldSeparation( previousFragIid, frag->iid, &fragsScaffoldSeparation, TRUE) != -1)
          {  
            fprintf( stderr, " fragsScaffoldSeparation: " F_COORD "\n", fragsScaffoldSeparation);
            fprintf( stderr, "   fragsLocaleSeparation: " F_COORD "\n", frag->localePos.bgn - previousFrag->localePos.end);
            sizesContigContainedGaps[ numContigContainedGaps++ ] = 
              abs( fragsScaffoldSeparation -
                   (frag->localePos.bgn - previousFrag->localePos.end));
            numContigContainedGapsPerLocale[ frag->locale ]++;					 
          }
          else
          {
            fragsMisoriented = TRUE;
          }
        }
        if ( previousUnitigId != unitig->id && fragsMisoriented == FALSE)
        {
          numPureContigContainedGaps++;
          fprintf( stderr, "***--- contig contained gap (pure) ---*** (" F_CID "," F_CID ")\n", beginFragIid, endFragIid);
        }
      }
      
      if ( fragsMisoriented == TRUE )
        fprintf( stderr, "***--- fragsMisoriented  ---*** (" F_CID "," F_CID ")\n", beginFragIid, endFragIid);
      
      if ( fragsMisoriented == FALSE )
      {
        if ( previousUnitigId == unitig->id && gapSizeInFrags == 1 )
        {
          fprintf( stderr, "***--- unitig contained gap ---*** (" F_CID "," F_CID ")\n", beginFragIid, endFragIid);
          if (previousLocale == frag->locale)
          {
            int fragsScaffoldSeparation;
            
            if (getScaffoldSeparation( previousFragIid, frag->iid, &fragsScaffoldSeparation, TRUE) != -1)
            {  
              fprintf( stderr, " fragsScaffoldSeparation: " F_COORD "\n", fragsScaffoldSeparation);
              fprintf( stderr, "   fragsLocaleSeparation: " F_COORD "\n", frag->localePos.bgn - previousFrag->localePos.end);
              sizesUnitigContainedGaps[ numUnitigContainedGaps++ ] = 
                abs( fragsScaffoldSeparation -
                     (frag->localePos.bgn - previousFrag->localePos.end));
              numUnitigContainedGapsPerLocale[ frag->locale ]++;
            }
          }
        }
        else if ( previousUnitigId == unitig->id )
        {
          fprintf( stderr, "***--- unitig contained gap (multiple frags) ---*** (" F_CID ", " F_CID ")\n", beginFragIid, endFragIid);
          if (previousLocale == frag->locale)
          {
            int fragsScaffoldSeparation;
            
            if (getScaffoldSeparation( previousFragIid, frag->iid, &fragsScaffoldSeparation, TRUE) != -1)
            {  
              fprintf( stderr, " fragsScaffoldSeparation: " F_COORD "\n", fragsScaffoldSeparation);
              fprintf( stderr, "   fragsLocaleSeparation: " F_COORD "\n", frag->localePos.bgn - previousFrag->localePos.end);
              sizesUnitigContainedGaps[ numUnitigContainedGaps++ ] = 
                abs( fragsScaffoldSeparation -
                     (frag->localePos.bgn - previousFrag->localePos.end));
              numUnitigContainedGapsPerLocale[ frag->locale ]++;
            }
          }
        }
        fprintf( stderr, "\n");
      }
    }
    prior = fragsNotInAss[ ifrag ];
  }
  
  fprintf( stderr, "\n\n");
  fprintf( stderr, "                  numBACFrags: %9ld\n", numBACFrags);
  fprintf( stderr, "             numBACFragsInAss: %9ld (%.2f %%)\n",
           numBACFragsInAss, 
           100.0 * numBACFragsInAss / numBACFrags);
  fprintf( stderr, "      numBACFragsInSurrogates: %9ld (%.2f %%)\n",
           numBACFragsInSurrogates, 
           100.0 * numBACFragsInSurrogates / numBACFrags);
  fprintf( stderr, "   numBACFragsNotInSurrogates: %9ld (%.2f %%)\n",
           numBACFragsNotInSurrogates,
           100.0 * numBACFragsNotInSurrogates / numBACFrags);
  fprintf( stderr, "       nBFIA + nBFIS + nBFNIS: %9ld (%.2f %%)\n", 
           numBACFragsInAss + numBACFragsInSurrogates + numBACFragsNotInSurrogates,
           100.0 * (numBACFragsInAss + numBACFragsInSurrogates + numBACFragsNotInSurrogates) / numBACFrags);
  fprintf( stderr, "avg Length unitigLengthsBFNIS: %9.2lf (%d unitigs)\n", 
           unitigLengthsBFNIS / numUnitigsBFNIS, numUnitigsBFNIS);
  fprintf( stderr, "\n");
  
  fprintf( stderr, "            number of locales: %5" F_CIDP "\n", maxLocale);
  fprintf( stderr, "      numMissingBACFragRanges: %5d\n",
           numMissingBACFragRanges);
  fprintf( stderr, "  flankingFragsInSameScaffold: %5d (same locale but not same contig)\n", flankingFragsInSameScaffold);
  fprintf( stderr, "      flankingContigsAdjacent: %5d\n",
           flankingContigsAdjacent);
  fprintf( stderr, "\n");
  
  fprintf( stderr, "       numUnitigContainedGaps: %7d\n",
           numUnitigContainedGaps);
  fprintf( stderr, "       numContigContainedGaps: %7d\n",
           numContigContainedGaps);
  fprintf( stderr, "   numPureContigContainedGaps: %7d\n",
           numPureContigContainedGaps);  
  fprintf( stderr, "\n");
  
  fprintf( stderr, "      \t             numBACFragsNot\t numUnitig    \t numContig    \t numMissingBAC\n");
  fprintf( stderr, "locale\t numBACFrags   InSurrogates\t ContainedGaps\t ContainedGaps\t FragRanges\n");
  for (i = 1; i <= maxLocale; i++)
  {
    fprintf( stderr, "%3d     %7d\t %7d  (%5.2f %%)\t %5d\t\t %5d\t\t %5d\n", 
             i, numFragsInLocale[i], missingFragsPerLocale[i],
             (100.0 * missingFragsPerLocale[i]) / numFragsInLocale[i],
             numUnitigContainedGapsPerLocale[i], numContigContainedGapsPerLocale[i],
             numMissingBACFragRangesPerLocale[i]);
  }
  fprintf( stderr, "\n");
  
  dumpCelagram( gapsFlanked, flankingFragsInSameScaffold, "gapsFlanked");
  dumpCelagram( sizesUnitigContainedGaps, numUnitigContainedGaps, "unitigContainedGaps");
  dumpCelagram( sizesContigContainedGaps, numContigContainedGaps, "contigContainedGaps");
  dumpCelagram( &numUnitigContainedGapsPerLocale[1], maxLocale, "unitigContainedGapsPerLocale");
  dumpCelagram( &numContigContainedGapsPerLocale[1], maxLocale, "contigContainedGapsPerLocale");
  
  
  localeFragsInUnitig = (int *) safe_malloc ( GetNumGraphNodes( ScaffoldGraph->CIGraph ) * sizeof( int ));
  if (localeFragsInUnitig == NULL)
  {
    fprintf( stderr, "could not safe_malloc space for localeFragsInUnitig!\n");
    assert(0);
  }
  for (i = 0; i < GetNumGraphNodes( ScaffoldGraph->CIGraph ); i++)
    localeFragsInUnitig[ i ] = 0;
  
  localeFragsInContig = (int *) safe_malloc ( GetNumGraphNodes( ScaffoldGraph->ContigGraph ) * sizeof( int ));
  if (localeFragsInContig == NULL)
  {
    fprintf( stderr, "could not safe_malloc space for localeFragsInContig!\n");
    assert(0);
  }
  for (i = 0; i < GetNumGraphNodes( ScaffoldGraph->ContigGraph ); i++)
    localeFragsInContig[ i ] = 0;
  
  localeFragsInScaffold = (int *) safe_malloc ( GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph ) * sizeof( int ));
  if (localeFragsInScaffold == NULL)
  {
    fprintf( stderr, "could not safe_malloc space for localeFragsInScaffold!\n");
    assert(0);
  }
  for (i = 0; i < GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph ); i++)
    localeFragsInScaffold[ i ] = 0;
  
  // these loops really suck, but since frags from a locale are not dense in fragIid land
  // we don't have a lot of choice
  fprintf( stderr, "locale\t numUnitigsInLocale\t numContigsInLocale\t numScaffoldsInLocale\t missingFragsRecount\n");
  for (i = 1; i <= maxLocale; i++)
  {	
    for (ifrag = 0; ifrag < GetNumVA_CIFragT( ScaffoldGraph->CIFrags ); ifrag++)
    {
      frag = GetCIFragT( ScaffoldGraph->CIFrags, ifrag);
      
      if (frag->locale == i)
      {
        CDS_CID_t scaffID = GetGraphNode(ScaffoldGraph->ContigGraph, frag->contigID)->scaffoldID;
        // we need to do something with surrogates to get the true contig/scaffold stats
        // NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->cid);
        
        
        {
          CDS_CID_t tempiid = 370253;
          NodeCGW_T *unitiggy = GetGraphNode(ScaffoldGraph->CIGraph, frag->cid);
          
          if (frag->iid == tempiid)
            fprintf( stderr, "frag " F_CID " shows as being in scaff " F_CID ", unitig: " F_CID " (numInstances: %d)\n", 
                     tempiid, scaffID, frag->cid, unitiggy->info.CI.numInstances);
        }
        
        
        
        
        if ( localeFragsInUnitig[ frag->cid ] != i)
        {
          numUnitigsInLocale++;
          localeFragsInUnitig[ frag->cid ] = i;
        }
        if ( localeFragsInContig[ frag->contigID ] != i)
        {
          numContigsInLocale++;
          localeFragsInContig[ frag->contigID ] = i;
        }
        if ( localeFragsInScaffold[ scaffID ] != i && scaffID != NULLINDEX)
        {
          numScaffoldsInLocale++;
          localeFragsInScaffold[ scaffID ] = i;
          fprintf( stderr, "locale %d claims to be in scaffold " F_CID " via frag iid " F_CID " (unitig " F_CID ")\n", 
                   i, scaffID, frag->iid, frag->cid);
        }
        
        // some temp stuff for double checking
        if (scaffID == NULLINDEX)
          missingFragsRecount++;
        {
          NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->cid);
          if ( unitig->info.CI.numInstances > 0 )
            missingFragsRecount--;	
        }
      }
    }
    fprintf( stderr, "%d\t %7d\t\t %7d\t\t %7d\t\t %d\n", i, 
             numUnitigsInLocale, numContigsInLocale, numScaffoldsInLocale, missingFragsRecount);
    numUnitigsInLocale = 0;
    numContigsInLocale = 0;
    numScaffoldsInLocale = 0;	
    missingFragsRecount = 0;
  }
  
  
  {
    int numRealScaffolds = 0;
    
    for ( i= 0; i < GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph); i++)
    {
      NodeCGW_T *scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, i);
      if(scaff->type == REAL_SCAFFOLD)
        numRealScaffolds++;
    }
    fprintf( stderr, "numRealScaffolds: %d\n", numRealScaffolds);
  }
  
  prepareBACFragInfo( numBACFrags, maxLocale);
  
  exit(0);
}

int	getScaffoldSeparation( CDS_CID_t frag1Iid,
                               CDS_CID_t frag2Iid,
                               int *fragScaffoldSeparation,
                               int verbose)
{
  int i;
  CIFragT *frag1, *frag2;
  int scaffold1, scaffold2;
  CDS_COORD_t frag1LeftEnd, frag1RightEnd;
  CDS_COORD_t frag2LeftEnd, frag2RightEnd;
  int frag2Orientation, frag1Orientation;
  InfoByIID *info;
  
  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, frag1Iid);
  assert(info->set);
  frag1 = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
  
  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, frag2Iid);
  assert(info->set);
  frag2 = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
  
  GetFragmentPositionInScaffold( frag1, &frag1LeftEnd, &frag1RightEnd, &frag1Orientation);
  GetFragmentPositionInScaffold( frag2, &frag2LeftEnd, &frag2RightEnd, &frag2Orientation);
  
  if (frag1Orientation != frag2Orientation)
  {
    if (verbose)
      fprintf( stderr, "frags " F_CID " and " F_CID " have different orientations!\n", frag1Iid, frag2Iid);
    return -1;
  }
  
  scaffold1 = GetGraphNode( ScaffoldGraph->ContigGraph, frag1->contigID)->scaffoldID;
  scaffold2 = GetGraphNode( ScaffoldGraph->ContigGraph, frag2->contigID)->scaffoldID;
  
  if (scaffold1 != scaffold2)
  {
    if (verbose)
      fprintf( stderr, "frags " F_CID " and " F_CID " are in different scaffolds!\n", frag1Iid, frag2Iid);
    assert( 0 );
  }
  
  if ( min( frag1LeftEnd, frag1RightEnd) < min( frag2LeftEnd, frag2RightEnd))          // frag1 left of frag2
    *fragScaffoldSeparation = min( frag2LeftEnd, frag2RightEnd) - max( frag1LeftEnd, frag1RightEnd);
  else if ( min( frag2LeftEnd, frag2RightEnd) < min( frag1LeftEnd, frag1RightEnd))     // frag2 left of frag1
    *fragScaffoldSeparation = min( frag1LeftEnd, frag1RightEnd) - max( frag2LeftEnd, frag2RightEnd);
  
  return 1;
}

void dumpCelagram( int *array, int numArrayItems, char *label)
{
  char fname[1024];
  int i;
  FILE *outFile;
  
  sprintf( fname, "%s.cgm", label);
  outFile = fopen( fname, "w");
  if (outFile == NULL)
  {
    fprintf( stderr, "can not open file named %s in dumpCelagram\n", fname);
    exit(1);
  }
  
  fprintf( outFile, "%s\n", label);
  for (i = 0; i < numArrayItems; i++)
    fprintf( outFile, "%d\n", array[i]);
  
  fclose( outFile );
}

int compareLocales(const void *e1,
                   const void *e2) 
{
  // sort by locale, then by locale position
  
  fragDataT *frag1 = (fragDataT *) e1;
  fragDataT *frag2 = (fragDataT *) e2;
  
  assert((frag1 != NULL) && (frag2 != NULL));
  
  if ( frag1->locale < frag2->locale)
    return -1;
  else if ( frag1->locale > frag2->locale)
    return 1;
  else // sort by pos
  {
    if ( frag1->localePos.bgn < frag2->localePos.bgn )
      return -1;
    else if ( frag1->localePos.bgn > frag2->localePos.bgn )
      return 1;
    else
    {
      fprintf( stderr, "frag1->iid: " F_CID ", frag1->localePos.bgn: " F_COORD "\n",
               frag1->iid, frag1->localePos.bgn);
      fprintf( stderr, "frag2->iid: " F_CID ", frag2->localePos.bgn: " F_COORD "\n",
               frag2->iid, frag2->localePos.bgn);
      assert(0);
    }
  }
  return(0);  // useless statement to make the compiler stop complaining
}

void prepareBACFragInfo( int numBACFrags, CDS_CID_t maxLocale)
{
  int ifrag, ilocale, BACFragCount = 0;
  fragDataT *allBACFragsData;
  CIFragT *frag;
  int fragsPerLocale[ maxLocale ], localeLength[ maxLocale ];
  rangeDataT *contigConsistentRanges;
  int numContigConsistentRanges;
  
  allBACFragsData = ( fragDataT *) safe_malloc( (GetNumVA_CIFragT( ScaffoldGraph->CIFrags ) + 1) * sizeof( fragDataT ));
  if (allBACFragsData == NULL)
  {
    fprintf( stderr, "failed to safe_malloc space for allBACFragsData!\n");
    assert(0);
  }
  
  for ( ilocale = 1; ilocale <= maxLocale; ilocale++)
    fragsPerLocale[ ilocale ] = 0;
  
  for (ifrag = 0; ifrag < GetNumVA_CIFragT( ScaffoldGraph->CIFrags ); ifrag++)
  {
    CIFragT *frag = GetCIFragT( ScaffoldGraph->CIFrags, ifrag);
    if (frag->locale >= 0)
    {
      gatherBACFragInfo( frag, &allBACFragsData[ BACFragCount++ ]);
      fragsPerLocale[ frag->locale ]++;
      if (frag->localePos.end > localeLength[ frag->locale ])
        localeLength[ frag->locale ] = frag->localePos.end;
    }
  }
  
  // sort by locale, then by position in locale
  qsort( allBACFragsData, numBACFrags, sizeof( fragDataT ), compareLocales);
  
  if (0)
    for ( ifrag = 0; ifrag < numBACFrags; ifrag++)
    {
      fprintf( stderr, "allBACFragsData[%7d]: %2" F_CIDP " %7" F_COORDP " %7" F_COORDP,
               ifrag, allBACFragsData[ifrag].locale,
               allBACFragsData[ifrag].localePos.bgn,
               allBACFragsData[ifrag].localePos.end);
      if (ifrag > 0 && 
          allBACFragsData[ifrag].localePos.bgn > allBACFragsData[ifrag - 1].localePos.end &&
          allBACFragsData[ifrag].locale == allBACFragsData[ifrag - 1].locale)
        fprintf( stderr, " ***** delta ");
      fprintf( stderr, "\n");
    }
  
  computeGoodnessMeasure( allBACFragsData, numBACFrags, fragsPerLocale, maxLocale);
  
  contigConsistentRanges = (rangeDataT *) safe_malloc( numBACFrags * sizeof( rangeDataT ));
  if ( contigConsistentRanges == NULL)
  {
    fprintf( stderr, "could not safe_malloc space (" F_SIZE_T " bytes) for contigConsistentRanges!\n",
             numBACFrags * sizeof( rangeDataT ));
    assert(0);
  }
  
  numContigConsistentRanges = computeConsistentRanges( allBACFragsData, numBACFrags, 
                                                       fragsPerLocale, maxLocale, 
                                                       contigConsistentRanges);
  
  dumpGraphFile( allBACFragsData, numBACFrags,
                 contigConsistentRanges, numContigConsistentRanges,
                 localeLength, maxLocale);
  
  compareBACsToScaffolds( allBACFragsData, numBACFrags, fragsPerLocale, maxLocale, 
                          contigConsistentRanges, numContigConsistentRanges);
  
  free ( contigConsistentRanges );
  free ( allBACFragsData );  // need to step through and free surrogates as well
}

void computeGoodnessMeasure( fragDataT *allBACFragsData,
                             int numBACFrags,
                             int *fragsPerLocale,
                             CDS_CID_t maxLocale)
{
  int ifrag, ilocale;
  int unitigCorrectOverlaps[ maxLocale ], contigCorrectOverlaps[ maxLocale ], scaffoldCorrectOverlaps[ maxLocale ];
  
  for( ilocale = 1; ilocale <= maxLocale; ilocale++)
    unitigCorrectOverlaps[ ilocale ] = contigCorrectOverlaps[ ilocale ] = scaffoldCorrectOverlaps[ ilocale ] = 0;
  
  for ( ifrag = 0; ifrag < numBACFrags - 1; ifrag++)
  {
    int localeOverlap = abs( allBACFragsData[ifrag].localePos.end - allBACFragsData[ifrag + 1].localePos.bgn);
    int unitigOverlap = abs( allBACFragsData[ifrag].unitigOffset3p - allBACFragsData[ifrag + 1].unitigOffset5p);
    int contigOverlap = abs( allBACFragsData[ifrag].contigOffset3p - allBACFragsData[ifrag + 1].contigOffset5p);
    int scaffoldOverlap = abs( allBACFragsData[ifrag].scaffoldOffset3p - allBACFragsData[ifrag + 1].scaffoldOffset5p);
    int unitigCorrect, contigCorrect, scaffoldCorrect;
    
    unitigCorrect = contigCorrect = scaffoldCorrect = FALSE;
    
    if ( abs( localeOverlap - unitigOverlap ) < UNITIG_CONSISTENCY_CONSTANT && 
         allBACFragsData[ifrag].unitigID == allBACFragsData[ifrag + 1].unitigID &&
         allBACFragsData[ifrag].locale == allBACFragsData[ifrag + 1].locale)
    { 
      unitigCorrectOverlaps[ allBACFragsData[ifrag].locale ]++;
      unitigCorrect = TRUE;
    }
    
    if ( abs( localeOverlap - contigOverlap ) < CONTIG_CONSISTENCY_CONSTANT && 
         allBACFragsData[ifrag].contigID == allBACFragsData[ifrag + 1].contigID &&
         allBACFragsData[ifrag].locale == allBACFragsData[ifrag + 1].locale)
    {
      contigCorrectOverlaps[ allBACFragsData[ifrag].locale ]++;
      contigCorrect = TRUE;
    }
    
    // since unscaffolded contigs all have scaffoldID = NULLINDEX have to dig a bit
    if ( allBACFragsData[ifrag].scaffoldID == allBACFragsData[ifrag + 1].scaffoldID)
    {
      // the case where they are in the same real scaffold
      if ( abs( localeOverlap - scaffoldOverlap ) < SCAFFOLD_CONSISTENCY_CONSTANT && 
           allBACFragsData[ifrag].scaffoldID != NULLINDEX &&
           allBACFragsData[ifrag].locale == allBACFragsData[ifrag + 1].locale)
      {
        scaffoldCorrectOverlaps[ allBACFragsData[ifrag].locale ]++;
        scaffoldCorrect = TRUE;
      }
      
      // the case where they are in the same contig but the scaffold has NULLINDEX
      if ( contigCorrect == TRUE &&
           allBACFragsData[ifrag].scaffoldID == NULLINDEX &&
           allBACFragsData[ifrag].locale == allBACFragsData[ifrag + 1].locale)
      {
        scaffoldCorrectOverlaps[ allBACFragsData[ifrag].locale ]++;
        scaffoldCorrect = TRUE;
      }
    }
    
    if ( !contigCorrect && scaffoldCorrect )
      fprintf( stderr, "overlap in scaffold but not in contig for frags " F_CID " and " F_CID " (locale: " F_CID ")\n",
               allBACFragsData[ifrag].iid, allBACFragsData[ifrag + 1].iid, allBACFragsData[ifrag + 1].locale);
    
  }
  
  fprintf( stderr,   "                                         unitig                contig             scaff\n");
  fprintf( stderr,   "locale     fragsPerLocale       CorrectOverlaps       CorrectOverlaps   CorrectOverlaps\n");
  for( ilocale = 1; ilocale <= maxLocale; ilocale++)
    fprintf( stderr, "%6d            %7d                 %5d                 %5d             %5d\n", 
             ilocale, fragsPerLocale[ ilocale ], 
             unitigCorrectOverlaps[ ilocale ],
             contigCorrectOverlaps[ ilocale ],
             scaffoldCorrectOverlaps[ ilocale ]);  
}

// this computes the consistent ranges on a contig basis
int computeConsistentRanges( fragDataT *allBACFragsData,
                             int numBACFrags, 
                             int *fragsPerLocale,
                             CDS_CID_t maxLocale, 
                             rangeDataT *contigConsistentRanges)
{
  int i, icnt, ifrag, ilocale, inRange = FALSE;
  int totalFrags = 0, numConsistentRanges = 0;
  int ifragStart;
  CDS_CID_t currentContigID = NULLINDEX;
  fragDataT *allBACFragsDataWithSurrogates;
  
  // first construct a version of allBACFragsData that has the surrogate appearances 
  // of frags on an equal footing with the real ones, but no parents of surrogates
  // count them
  for ( ifrag = 0; ifrag < numBACFrags; ifrag++)
  {
    if (allBACFragsData[ ifrag ].surrogateCount > 0)
      totalFrags += allBACFragsData[ ifrag ].surrogateCount;  // count each surrogate
    else
      totalFrags++;  // add original
    
    if (totalFrags < ifrag)
      fprintf( stderr, "WTF?\n");
  }
  
  // make space for them
  allBACFragsDataWithSurrogates = (fragDataT *) safe_malloc( totalFrags * sizeof( fragDataT ));
  if ( allBACFragsDataWithSurrogates == NULL)
    assert(0);
  
  // transfer them
  icnt = 0;
  for ( ifrag = 0; ifrag < numBACFrags; ifrag++)
  {
    if ( allBACFragsData[ ifrag ].surrogateCount == 0)  // frag exists only in a unitig, no surrogates
      allBACFragsDataWithSurrogates[ icnt++ ] = allBACFragsData[ ifrag ];
    else  // frag exists in at least one surrogate
      for (i = 0; i < allBACFragsData[ ifrag ].surrogateCount; i++)
        allBACFragsDataWithSurrogates[ icnt++ ] = allBACFragsData[ ifrag ].surrogates[i];
    
    if (allBACFragsDataWithSurrogates[ icnt - 1].scaffoldOffset3p < 0)
      fprintf( stderr, "just copied some crap!\n");
  }
  
  // now sort by contig, secondarily by contigOffset5p
  qsort( allBACFragsDataWithSurrogates, totalFrags, sizeof( fragDataT ), compContigIDs);
  
  // now step through contigs, finding the consistent ranges
  inRange = FALSE;
  for ( ifrag = 0; ifrag < totalFrags - 1; ifrag++)
  {
    int localeOverlap, contigOverlap;
    
    if ( allBACFragsDataWithSurrogates[ifrag].iid == 370155 ||
         allBACFragsDataWithSurrogates[ifrag].iid == 370150)
      fprintf( stderr, "at problem ifrag\n");
    
    // now compare the positions of ifragData and ifragPlusOneData 
    localeOverlap = abs( allBACFragsDataWithSurrogates[ifrag].localePos.end - 
                         allBACFragsDataWithSurrogates[ifrag + 1].localePos.bgn);
    contigOverlap = abs( allBACFragsDataWithSurrogates[ifrag].contigOffset3p - 
                         allBACFragsDataWithSurrogates[ifrag + 1].contigOffset5p );
    
    if ( inRange == FALSE )
    {
      fprintf( stderr, "begin: " F_CID " (ifrag: %d)\n", allBACFragsDataWithSurrogates[ifrag].iid, ifrag);
      inRange = TRUE;
      currentContigID = allBACFragsDataWithSurrogates[ifrag].contigID;
      contigConsistentRanges[ numConsistentRanges ].locale = allBACFragsDataWithSurrogates[ifrag].locale;
      contigConsistentRanges[ numConsistentRanges ].scaffoldID = allBACFragsDataWithSurrogates[ifrag].scaffoldID;
      contigConsistentRanges[ numConsistentRanges ].contigID = allBACFragsDataWithSurrogates[ifrag].contigID;
      contigConsistentRanges[ numConsistentRanges ].iidStart = allBACFragsDataWithSurrogates[ifrag].iid;
      // contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = 
      // allBACFragsDataWithSurrogates[ifrag].scaffoldOffset5p;	  
      // contigConsistentRanges[ numConsistentRanges ].localePosStart = allBACFragsDataWithSurrogates[ifrag].localePos;
      
      ifragStart = ifrag;
    }
    
    if ( allBACFragsDataWithSurrogates[ifrag].iid >= 486916 && 
         allBACFragsDataWithSurrogates[ifrag].iid <= 486954)
    {
      fprintf( stderr, "anfjdanf allBACFragsDataWithSurrogates[ %d ].contigID (iid: " F_CID "): " F_CID "\n",
               ifrag, allBACFragsDataWithSurrogates[ifrag].iid, allBACFragsDataWithSurrogates[ifrag].contigID);
    }
    
    if ( abs( localeOverlap - contigOverlap ) < CONTIG_CONSISTENCY_CONSTANT && // overlap properly
         abs( allBACFragsDataWithSurrogates[ifrag].iid - 
              allBACFragsDataWithSurrogates[ifrag + 1].iid == 1) &&       
         allBACFragsDataWithSurrogates[ifrag].contigID == allBACFragsDataWithSurrogates[ifrag + 1].contigID &&
         allBACFragsDataWithSurrogates[ifrag].locale == allBACFragsDataWithSurrogates[ifrag + 1].locale)
    {
      // fprintf( stderr, F_CID " <-> " F_CID " are consistent\n", allBACFragsDataWithSurrogates[ifrag].iid, 
      // allBACFragsDataWithSurrogates[ifrag + 1].iid);
    }
    else
    {
      if (inRange)
      {
        NodeCGW_T *contig = GetGraphNode( ScaffoldGraph->ContigGraph, allBACFragsDataWithSurrogates[ifrag].contigID);
        int contigOrientation = ( contig->offsetAEnd.mean < contig->offsetBEnd.mean ) ? 0 : 1;
        
        if ( allBACFragsDataWithSurrogates[ifrag].iid == 413190)
          fprintf( stderr, "\n");
        
        fprintf( stderr, "end: " F_CID " (ifrag: %d)\n", allBACFragsDataWithSurrogates[ifrag].iid, ifrag);
        inRange = FALSE;
        currentContigID = -1;	
        contigConsistentRanges[ numConsistentRanges ].iidEnd = allBACFragsDataWithSurrogates[ifrag].iid;
        
        contigConsistentRanges[ numConsistentRanges ].localePosStart = 
          allBACFragsDataWithSurrogates[ ifragStart ].localePos;
        contigConsistentRanges[ numConsistentRanges ].localePosEnd = 
          allBACFragsDataWithSurrogates[ ifrag ].localePos;
        
        if ( contigConsistentRanges[ numConsistentRanges ].iidStart < 
             contigConsistentRanges[ numConsistentRanges ].iidEnd )       // ascending in frag iid in the contig
        {
          if (contigOrientation == 0)
          {
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = 
              allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset5p;	  
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMax = 
              allBACFragsDataWithSurrogates[ ifrag ].scaffoldOffset3p;	  
          }
          else // contigOrientation == 1
          {
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = 
              allBACFragsDataWithSurrogates[ ifrag ].scaffoldOffset3p;	  
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMax = 
              allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset5p;	  
          }
        }
        else if ( contigConsistentRanges[ numConsistentRanges ].iidStart > 
                  contigConsistentRanges[ numConsistentRanges ].iidEnd )  // descending in iid in the contig
        {
          if (contigOrientation == 0)
          {
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = 
              allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset3p;	  
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMax = 
              allBACFragsDataWithSurrogates[ ifrag ].scaffoldOffset5p;	  
          }
          else // contigOrientation == 1
          {
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = 
              allBACFragsDataWithSurrogates[ ifrag ].scaffoldOffset5p;	  
            contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMax = 
              allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset3p;	  
          }
        }		  
        else  // single frag in range
        {
          contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = 
            min( allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset3p,
                 allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset5p);
          contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMax = 
            max( allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset3p,
                 allBACFragsDataWithSurrogates[ ifragStart ].scaffoldOffset5p);
        }
        
#if 0
        if ( contigConsistentRanges[ numConsistentRanges ].iidEnd - contigConsistentRanges[ numConsistentRanges ].iidStart 
             > 10 
             &&
             ( contigConsistentRanges[ numConsistentRanges ].iidStart >
               contigConsistentRanges[ numConsistentRanges ].iidEnd ) )
          fprintf( stderr, "long stretch!\n");
#endif		  
        
        if ( contigConsistentRanges[ numConsistentRanges ].iidEnd == 413190)
          fprintf( stderr, "at contig range ending with iid %d\n", 413190);
        
        numConsistentRanges++;
      }
    }
  }
  
  if (1)
  {
    fprintf( stderr, "range    iidStart      iidEnd    locale     scaffoldID    numFrags\n");
    for (i = 0; i < numConsistentRanges; i++)
    {
      fprintf( stderr, "%5d    %8d      %6d    %6d     %10d    %8d\n", i, contigConsistentRanges[i].iidStart, 
               contigConsistentRanges[i].iidEnd, contigConsistentRanges[i].locale, 
               contigConsistentRanges[i].scaffoldID, 
               abs( contigConsistentRanges[i].iidEnd - contigConsistentRanges[i].iidStart) + 1);
      
      if ( contigConsistentRanges[i].scaffoldOffsetMin > contigConsistentRanges[i].scaffoldOffsetMax)
      {
        fprintf( stderr, "contigConsistentRanges[ %d ].scaffoldOffsetMin: %f\n",
                 i, contigConsistentRanges[i].scaffoldOffsetMin);
        fprintf( stderr, "contigConsistentRanges[ %d ].scaffoldOffsetMax: %f\n",
                 i, contigConsistentRanges[i].scaffoldOffsetMax);
        assert(0);
      }
    }
  }
  free ( allBACFragsDataWithSurrogates );
  return( numConsistentRanges );
}

// this computes the consistent ranges on a contig basis
int computeConsistentRanges_old( fragDataT *allBACFragsData,
                                 int numBACFrags, 
                                 int *fragsPerLocale,
                                 CDS_CID_t maxLocale, 
                                 rangeDataT *contigConsistentRanges)
{
  int i, ifrag, ilocale, inRange = FALSE;
  int numConsistentRanges = 0;
  fragDataT *closestFragData;
  CDS_CID_t currentContigID = NULLINDEX;
  fragDataT *ifragData, *ifragPlusOneData;
  
  for ( ifrag = 0; ifrag < numBACFrags - 1; ifrag++)
  {
    int localeOverlap = abs( allBACFragsData[ifrag].localePos.end - allBACFragsData[ifrag + 1].localePos.bgn);
    int contigOverlap;
    int contigCorrect;
    
    // initialize pointers to the base fragDataT for both ifrag and ifrag+1
    ifragData = &allBACFragsData[ifrag];
    ifragPlusOneData = &allBACFragsData[ifrag + 1];
    // contigOverlap = abs( allBACFragsData[ifrag].contigOffset3p - allBACFragsData[ifrag + 1].contigOffset5p);
    
    // if we are currently in a contig and ifrag has surrogates find the surrogate that's in this contig
    if ( allBACFragsData[ifrag].surrogateCount > 0 && currentContigID != -1 )
    {
      for ( i = 0; i < allBACFragsData[ifrag].surrogateCount; i++)
      {
        if ( allBACFragsData[ifrag].surrogates[i].contigID == currentContigID )
        {		  
          ifragData = &(allBACFragsData[ifrag].surrogates[i]);
          break;
        }
      }
    }
    
    // if we are currently in a contig and ifrag+1 has surrogates find the surrogate that's in this contig
    if ( allBACFragsData[ifrag + 1].surrogateCount > 0 && currentContigID != -1 )
    {
      for ( i = 0; i < allBACFragsData[ifrag + 1].surrogateCount; i++)
      {
        if ( allBACFragsData[ifrag + 1].surrogates[i].contigID == currentContigID )
        {		  
          ifragPlusOneData = &(allBACFragsData[ifrag + 1].surrogates[i]);
          break;
        }
      }
    }
    
    // now compare the positions of ifragData and ifragPlusOneData 
    contigOverlap = abs( ifragData->contigOffset3p - ifragPlusOneData->contigOffset5p );
    
    if ( allBACFragsData[ifrag].iid == 412443 ||
         allBACFragsData[ifrag].iid == 412446)
      fprintf( stderr, "frag " F_CID " has contigOffset3p %d in contig " F_CID " contigOverlap %d localeOverlap %d\n",
               allBACFragsData[ifrag].iid, allBACFragsData[ifrag].contigOffset3p, allBACFragsData[ifrag].contigID,
               contigOverlap, localeOverlap);
    
    contigCorrect = FALSE;
    
    if ( inRange == FALSE && allBACFragsData[ifrag].surrogateCount == 0)
    {
      fprintf( stderr, "begin: " F_CID "\n", allBACFragsData[ifrag].iid);
      inRange = TRUE;
      currentContigID = allBACFragsData[ifrag].contigID;
      contigConsistentRanges[ numConsistentRanges ].locale = allBACFragsData[ifrag].locale;
      contigConsistentRanges[ numConsistentRanges ].scaffoldID = allBACFragsData[ifrag].scaffoldID;
      contigConsistentRanges[ numConsistentRanges ].contigID = allBACFragsData[ifrag].contigID;
      contigConsistentRanges[ numConsistentRanges ].iidStart = allBACFragsData[ifrag].iid;
      contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMin = allBACFragsData[ifrag].scaffoldOffset5p;	  
    }
    
    if ( abs( localeOverlap - contigOverlap ) < CONTIG_CONSISTENCY_CONSTANT && 
         ifragData->contigID == ifragPlusOneData->contigID &&
         ifragData->locale == ifragPlusOneData->locale)
    {
      // fprintf( stderr, F_CID " <-> " F_CID " are consistent\n", allBACFragsData[ifrag].iid, allBACFragsData[ifrag + 1].iid);
    }
    else
    {
      if (inRange)
      {
        fprintf( stderr, "end: " F_CID "\n", allBACFragsData[ifrag].iid);
        inRange = FALSE;
        currentContigID = NULLINDEX;	
        contigConsistentRanges[ numConsistentRanges ].iidEnd = ifragData->iid;
        contigConsistentRanges[ numConsistentRanges ].scaffoldOffsetMax = ifragData->scaffoldOffset3p;	  
        numConsistentRanges++;
      }
    }
  }
  
  if (1)
  {
    fprintf( stderr, "range    iidStart      iidEnd    locale     scaffoldID    numFrags\n");
    for (i = 0; i < numConsistentRanges; i++)
    {
      fprintf( stderr, "%5d    %8d      %6d    %6d     %10d    %8d\n", i, contigConsistentRanges[i].iidStart, 
               contigConsistentRanges[i].iidEnd, contigConsistentRanges[i].locale, 
               contigConsistentRanges[i].scaffoldID, 
               contigConsistentRanges[i].iidEnd - contigConsistentRanges[i].iidStart + 1);
    }
  }
  return( numConsistentRanges );
}



// what we have to do here is gather all the ranges in a scaffold
// then sort them by position
// then we can call computeBACScaffoldConsistency

void compareBACsToScaffolds( fragDataT *allBACFragsData,
                             int numBACFrags,
                             int *fragsPerLocale,
                             CDS_CID_t maxLocale, 
                             rangeDataT *contigConsistentRanges,
                             int numContigConsistentRanges)
{
  int i, ifrag, ilocale, irange, iscaff;
  int numScaffolds = GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph );
  scaffoldRangeT *BACScaffoldRange = (scaffoldRangeT *) safe_malloc ( numScaffolds * sizeof( scaffoldRangeT ));
  
  // sort contigConsistentRanges by scaffold, then by location in the scaffold
  qsort ( contigConsistentRanges, numContigConsistentRanges, sizeof( rangeDataT ), &compConsistentRangesByScaffold);
  
  // for each locale find the scaffolds that are in that locale
  for ( ilocale = 0; ilocale <= maxLocale; ilocale++)
  {
    // clear mapping of ranges to scaffolds for current locale
    for ( i = 0; i < numScaffolds; i++)
      BACScaffoldRange[i].firstRange = -1;
    
    // now step through ranges, seeing which scaffolds map to current locale 
    // when we're done, BACScaffoldRange[iscaffold] will hold the first and last contig consistent range
    // for scaffold "iscaffold" as determined by position in the scaffold
    for ( irange = 0; irange < numContigConsistentRanges; irange++)
    {
      if ( contigConsistentRanges[ irange ].scaffoldID != NULLINDEX &&
           contigConsistentRanges[ irange ].locale == ilocale)
      {
        // store the starting range in BACScaffoldRange.firstRange
        if ( BACScaffoldRange[ contigConsistentRanges[ irange ].scaffoldID ].firstRange == -1)
        {
          BACScaffoldRange[ contigConsistentRanges[ irange ].scaffoldID ].scaffoldID = 
            contigConsistentRanges[ irange ].scaffoldID;
          BACScaffoldRange[ contigConsistentRanges[ irange ].scaffoldID ].firstRange = irange;
          BACScaffoldRange[ contigConsistentRanges[ irange ].scaffoldID ].lastRange = irange;
        }
        else  // we also want the ending range, which is last irange that hits
        {
          BACScaffoldRange[ contigConsistentRanges[ irange ].scaffoldID ].lastRange = irange;
        }
      }
    }
    
    // now we have all the ranges for all scaffolds in this locale
    for ( iscaff = 0; iscaff < numScaffolds; iscaff++)
    {
      if ( BACScaffoldRange[ iscaff ].firstRange != -1 )
        computeBACScaffoldConsistency( iscaff, ilocale, allBACFragsData, numBACFrags,
                                       contigConsistentRanges, numContigConsistentRanges,
                                       BACScaffoldRange[ iscaff ].firstRange,
                                       BACScaffoldRange[ iscaff ].lastRange);
    }
    
    findContainedScaffolds( contigConsistentRanges, numContigConsistentRanges, BACScaffoldRange, ilocale );
  }
}

void findContainedScaffolds( rangeDataT *contigConsistentRanges,
                             int numContigConsistentRanges,
                             scaffoldRangeT *BACScaffoldRange,
                             CDS_CID_t locale )
{
  int i, j, numBACScaffolds = 0;
  int numScaffolds = GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph );
  scaffoldRangeT compressedScaffolds[ numScaffolds ];
  scaffoldIidRangeT compressedScaffoldsIids[ numScaffolds ];
  
  // first thin out the scaffolds not in this BAC
  for (i = 0; i < numScaffolds; i++)
  {
    if ( BACScaffoldRange[i].firstRange != -1)
      compressedScaffolds[ numBACScaffolds++ ] = BACScaffoldRange[i];
  }
  
  // now step through the compressed scaffolds, finding the first and last iids for each scaffold
  for ( i = 0; i < numBACScaffolds; i++)
  {
    // initialize to a valid value
    compressedScaffoldsIids[i].scaffoldID = compressedScaffolds[i].scaffoldID;
    compressedScaffoldsIids[i].iidStart = CDS_CID_MAX;
    compressedScaffoldsIids[i].iidEnd = NULLINDEX;
    
    fprintf( stderr, "compressedScaffolds[ %d ].firstRange: %d, compressedScaffolds[i].lastRange: %d\n",
             i, compressedScaffolds[i].firstRange, compressedScaffolds[i].lastRange);
    
    for ( j = compressedScaffolds[i].firstRange; j <= compressedScaffolds[i].lastRange; j++)
    {
      if (contigConsistentRanges[j].locale != locale)
        continue;
      
      // contigConsistentRanges[].iidStart is not necessarily less than contigConsistentRanges[].iidEnd
      // so have to check them both
      if ( contigConsistentRanges[j].iidStart < compressedScaffoldsIids[i].iidStart )
        compressedScaffoldsIids[i].iidStart = contigConsistentRanges[j].iidStart;
      if ( contigConsistentRanges[j].iidEnd < compressedScaffoldsIids[i].iidStart )
        compressedScaffoldsIids[i].iidStart = contigConsistentRanges[j].iidEnd;
      
      if ( contigConsistentRanges[j].iidStart > compressedScaffoldsIids[i].iidEnd )
        compressedScaffoldsIids[i].iidEnd = contigConsistentRanges[j].iidStart;
      if ( contigConsistentRanges[j].iidEnd > compressedScaffoldsIids[i].iidEnd )
        compressedScaffoldsIids[i].iidEnd = contigConsistentRanges[j].iidEnd;
    }
  }
  
  // now sort compressedScaffoldsIids by iid
  qsort ( compressedScaffoldsIids, numBACScaffolds, sizeof( scaffoldIidRangeT ), &compScaffoldIids);
  
  for ( i = 0; i < numBACScaffolds; i++)
  {
    fprintf( stderr, "compressedScaffoldsIids[ %d ].iidStart: %d, iidEnd: %d\n",
             i, compressedScaffoldsIids[i].iidStart, compressedScaffoldsIids[i].iidEnd);
    
    for ( j = i + 1; j < numBACScaffolds; j++)
    {
      if ( compressedScaffoldsIids[i].iidStart < compressedScaffoldsIids[j].iidStart &&
           compressedScaffoldsIids[i].iidEnd > compressedScaffoldsIids[j].iidEnd)
      {
        fprintf( stderr, 
                 "scaffold " F_CID " ( iids (%d, %d)) contains scaffold " F_CID " ( iids (%d, %d)), locale " F_CID "\n",
                 compressedScaffoldsIids[i].scaffoldID,
                 compressedScaffoldsIids[i].iidStart,
                 compressedScaffoldsIids[i].iidEnd,
                 compressedScaffoldsIids[j].scaffoldID,
                 compressedScaffoldsIids[j].iidStart,
                 compressedScaffoldsIids[j].iidEnd,
                 locale);
      }
    }
  }
}

int computeBACScaffoldConsistency( CDS_CID_t scaffoldID,
                                   CDS_CID_t locale, 
                                   fragDataT *allBACFragsData,
                                   int numBACFrags,
                                   rangeDataT *contigConsistentRanges,
                                   int numContigConsistentRanges,
                                   int firstRange,
                                   int lastRange)
{
  int irange, previousRange, tolerance = 25;
  CIFragT *previousFrag, *currFrag;
  int numRanges = 1, numInconsistentRanges = 0;
  double previousOffsetMin;
  InfoByIID *info;
  
  previousOffsetMin = contigConsistentRanges[ firstRange ].scaffoldOffsetMin;
  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, contigConsistentRanges[ firstRange ].iidStart);
  assert(info->set);
  previousFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
  previousRange = firstRange;
  
  for ( irange = firstRange + 1; irange <= lastRange; irange++)
  {
    int scaffPosDiff, localePosDiff;
    
    if ( contigConsistentRanges[ irange ].locale == locale && 
         contigConsistentRanges[ irange ].scaffoldID == scaffoldID && 
         contigConsistentRanges[ irange ].scaffoldID != NULLINDEX )
    {
      numRanges++;
      info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, contigConsistentRanges[ irange ].iidStart);
      assert(info->set);
      currFrag = GetCIFragT( ScaffoldGraph->CIFrags, info->fragIndex);
      
      scaffPosDiff = abs( previousOffsetMin - contigConsistentRanges[ irange ].scaffoldOffsetMin);
      localePosDiff = abs( currFrag->localePos.bgn - previousFrag->localePos.bgn );
      
      if ( abs( scaffPosDiff - localePosDiff) > tolerance)
        numInconsistentRanges++;
      
      previousFrag = currFrag;
      previousOffsetMin = contigConsistentRanges[ irange ].scaffoldOffsetMin;
      previousRange = irange;
    }
  }
  
  fprintf( stderr, "locale: %3" F_CIDP ", scaffoldID: %6" F_CIDP ", numRanges: %6d, numInconsistentRanges: %6d tolerance: %3d\n",
           locale, scaffoldID, numRanges, numInconsistentRanges, tolerance);
  
  return 1;
}


int gatherBACFragInfo( CIFragT *frag, fragDataT *allFragsData)
{
  NodeCGW_T *unitig, *contig, *scaffold;
  int i, numUnitigSurrogates;
  CDS_COORD_t fragLeftEnd, fragRightEnd;
  int fragOrientation;
  
  // first find out what the containing unitig is like, it determines the number of fragDataT's for this frag
  unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);
  assert (unitig != 0);  
  allFragsData->surrogateCount = unitig->info.CI.numInstances;
  
  contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
  assert (contig != 0);  
  
  // get info that's present for all frags regardless of surrogate involvement
  allFragsData->unitigID = frag->cid;
  allFragsData->unitigOffset5p = frag->offset5p.mean;
  allFragsData->unitigOffset3p = frag->offset3p.mean;	
  
  allFragsData->contigID = frag->contigID;
  allFragsData->contigOffset5p = frag->contigOffset5p.mean;
  allFragsData->contigOffset3p = frag->contigOffset3p.mean;
  
  allFragsData->iid = frag->iid;
  allFragsData->scaffoldID = contig->scaffoldID;
  allFragsData->locale = frag->locale;
  allFragsData->localePos = frag->localePos;
  
  if ( allFragsData->surrogateCount == 0 ) // the non-surrogate case
  {
    if ( allFragsData->scaffoldID != NULLINDEX)
    {
      GetFragmentPositionInScaffold( frag, &fragLeftEnd, &fragRightEnd, &fragOrientation);
      allFragsData->scaffoldID = contig->scaffoldID;
      if ( fragOrientation == 0 )
      {
        allFragsData->scaffoldOffset5p = fragLeftEnd;
        allFragsData->scaffoldOffset3p = fragRightEnd;
      }
      else
      {
        allFragsData->scaffoldOffset5p = fragRightEnd;
        allFragsData->scaffoldOffset3p = fragLeftEnd;
      }
    }
  }
  else // we have to step through surrogates
  {
    NodeCGW_T *surrogate;
    int instances[10];
    
    // put the surrogate info in an array
    allFragsData->surrogates = 
      (fragDataT *) safe_malloc ( allFragsData->surrogateCount * sizeof( fragDataT ));
    if (allFragsData->surrogates == NULL)
      assert(0);
    
    if (allFragsData->surrogateCount > 10)  // need to make instances[] bigger
      assert(0);
    
    // whether surrogateCount is greater than or less than 2 put it into instances array
    if (allFragsData->surrogateCount < 3)
    {
      instances[0] = unitig->info.CI.instances.in_line.instance1;
      if (allFragsData->surrogateCount == 2)
        instances[1] = unitig->info.CI.instances.in_line.instance2;
    }
    else
    {
      for ( i = 0; i < allFragsData->surrogateCount; i++)
        instances[i] = *GetVA_int( unitig->info.CI.instances.va, i);
    }
    
    // now that we have all the instance ids gathered, step through them collecting info
    for ( i = 0; i < allFragsData->surrogateCount; i++)
    {
      CDS_COORD_t surrogateLeftEnd, surrogateRightEnd;
      CDS_COORD_t contigOfSurrogateLeftEnd, contigOfSurrogateRightEnd;
      CDS_COORD_t fragLeftEnd, fragRightEnd;
      int surrogateScaffoldOrientation,
        contigOfSurrogateScaffoldOrientation,
        fragScaffoldOrientation;
      
      fragDataT *currSurrogateData = &(allFragsData->surrogates[i]);
      NodeCGW_T *contigOfSurrogate;
      
      surrogate = GetGraphNode( ScaffoldGraph->CIGraph, instances[i]);
      
      // frag's position in surrograte is same as in original untitig
      currSurrogateData->iid = frag->iid;
      currSurrogateData->locale = frag->locale;
      currSurrogateData->localePos = frag->localePos;
      currSurrogateData->unitigID = instances[i];
      currSurrogateData->unitigOffset5p = frag->offset5p.mean;
      currSurrogateData->unitigOffset3p = frag->offset3p.mean;
      
      currSurrogateData->contigID = surrogate->info.CI.contigID;
      contigOfSurrogate = GetGraphNode( ScaffoldGraph->ContigGraph, currSurrogateData->contigID);
      
#if 0
      GetContigPositionInScaffold( contigOfSurrogate, &contigOfSurrogateLeftEnd, 
                                   &contigOfSurrogateRightEnd, &contigOfSurrogateScaffoldOrientation);
      
      if ( contigOfSurrogate->offsetAEnd.mean < contigOfSurrogate->offsetBEnd.mean)  // contigOfSurrogate is A_B in scaffold
      {
        // determine surrogate orientation in contig
        if ( surrogate->offsetAEnd.mean < surrogate->offsetBEnd.mean )  // surrogate is A_B in contig
        {
          currSurrogateData->contigOffset5p = contigOfSurrogateLeftEnd + surrogate->offsetAEnd.mean + frag->offset5p.mean;
          currSurrogateData->contigOffset3p = contigOfSurrogateLeftEnd + surrogate->offsetAEnd.mean + frag->offset3p.mean;
        }
        else  // surrogate is flipped in contig
        {
          currSurrogateData->contigOffset5p = contigOfSurrogateLeftEnd + surrogate->offsetAEnd.mean - frag->offset5p.mean;
          currSurrogateData->contigOffset3p = contigOfSurrogateLeftEnd + surrogate->offsetAEnd.mean - frag->offset3p.mean;
        }
      }
      else  // contigOfSurrogate is B_A in its scaffold
      {
        // determine surrogate orientation in contig
        if ( surrogate->offsetAEnd.mean < surrogate->offsetBEnd.mean )  // surrogate is A_B in contig
        {
          currSurrogateData->contigOffset5p = contigOfSurrogateRightEnd - surrogate->offsetAEnd.mean - frag->offset5p.mean;
          currSurrogateData->contigOffset3p = contigOfSurrogateRightEnd - surrogate->offsetAEnd.mean - frag->offset3p.mean;
        }
        else  // surrogate is flipped in contig
        {
          currSurrogateData->contigOffset5p = contigOfSurrogateRightEnd - surrogate->offsetAEnd.mean + frag->offset5p.mean;
          currSurrogateData->contigOffset3p = contigOfSurrogateRightEnd - surrogate->offsetAEnd.mean + frag->offset3p.mean;
        }
      }
#endif
      
      // determine surrogate orientation in contig
      if ( surrogate->offsetAEnd.mean < surrogate->offsetBEnd.mean )  // surrogate is A_B in contig
      {
        currSurrogateData->contigOffset5p = surrogate->offsetAEnd.mean + frag->offset5p.mean;
        currSurrogateData->contigOffset3p = surrogate->offsetAEnd.mean + frag->offset3p.mean;
      }
      else  // surrogate is flipped in contig
      {
        currSurrogateData->contigOffset5p = surrogate->offsetAEnd.mean - frag->offset5p.mean;
        currSurrogateData->contigOffset3p = surrogate->offsetAEnd.mean - frag->offset3p.mean;
      }
      
      currSurrogateData->scaffoldID = contigOfSurrogate->scaffoldID;
      
      // now we need to compute frag's position in scaffold (via surrogate)
      if (frag->iid == 486952)
        GetCIPositionInScaffold( surrogate, &surrogateLeftEnd, &surrogateRightEnd, &surrogateScaffoldOrientation);
      else
        GetCIPositionInScaffold( surrogate, &surrogateLeftEnd, &surrogateRightEnd, &surrogateScaffoldOrientation);
      if (frag->iid == 486952)
        GetFragmentPositionInScaffoldFromCI( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation,
                                             surrogateLeftEnd, surrogateRightEnd, surrogateScaffoldOrientation);
      else
        GetFragmentPositionInScaffoldFromCI( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation,
                                             surrogateLeftEnd, surrogateRightEnd, surrogateScaffoldOrientation);
      
      if (fragScaffoldOrientation == 0)
      {
        currSurrogateData->scaffoldOffset5p = fragLeftEnd;
        currSurrogateData->scaffoldOffset3p = fragRightEnd;
      }	
      else
      {
        currSurrogateData->scaffoldOffset5p = fragRightEnd;
        currSurrogateData->scaffoldOffset3p = fragLeftEnd;
      }
    }
  }
  return 0;
}

void dumpGraphFile( fragDataT *allBACFragsData,
                    int numBACFrags,
                    rangeDataT *contigConsistentRanges,
                    int numContigConsistentRanges,
                    CDS_COORD_t *localeLength,
                    CDS_CID_t maxLocale)
{
  int i, irange;
  FILE *vizfile;
  
  vizfile = fopen("cam/scaff_BAC_matches", "w+");
  if (vizfile == NULL)
    assert(0);
  
  // first we dump scaffold lengths
  fprintf( vizfile, "{DATA pieces_x\n");
  for ( i = 0; i < GetNumGraphNodes( ScaffoldGraph->ScaffoldGraph ); i++)
  {
    NodeCGW_T *scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, i);
    
    if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD))
      continue;
    
    fprintf( vizfile, F_CID " 0: 1 %d\n", scaff->id, (int) scaff->bpLength.mean);
  }
  fprintf( vizfile, "}\n");
  
  // then we dump BAC lengths - do we have these somewhere?
  fprintf( vizfile, "{DATA pieces_y\n");
  for (i = 1; i <= maxLocale; i++)
    fprintf( vizfile, "%d 0: 1 " F_COORD "\n", i, localeLength[i]);	
  fprintf( vizfile, "}\n");
  
  // then we dump matches in the form scaff_start, bac_start, scaff_end, bac_end
  fprintf( vizfile, "{DATA matches_x_y\n");  
  for ( irange = 0; irange < numContigConsistentRanges; irange++)
  {
    if (contigConsistentRanges[ irange ].scaffoldID != NULLINDEX)
    {
      if ( contigConsistentRanges[ irange ].scaffoldOffsetMin < contigConsistentRanges[ irange ].scaffoldOffsetMax)
      {
        fprintf( vizfile, F_CID "_" F_CID " 0: " F_COORD " " F_COORD " " F_COORD " " F_COORD "\n",
                 contigConsistentRanges[ irange ].scaffoldID,
                 contigConsistentRanges[ irange ].locale,
                 (int) contigConsistentRanges[ irange ].scaffoldOffsetMin,
                 contigConsistentRanges[ irange ].localePosStart.bgn, 
                 (int) contigConsistentRanges[ irange ].scaffoldOffsetMax,
                 contigConsistentRanges[ irange ].localePosEnd.end);
      }
      else
        fprintf( vizfile, F_CID "_" F_CID " 0: " F_COORD " " F_COORD " " F_COORD " " F_COORD "\n",
                 contigConsistentRanges[ irange ].scaffoldID,
                 contigConsistentRanges[ irange ].locale,
                 (int) contigConsistentRanges[ irange ].scaffoldOffsetMax,
                 contigConsistentRanges[ irange ].localePosStart.bgn,
                 (int) contigConsistentRanges[ irange ].scaffoldOffsetMin,
                 contigConsistentRanges[ irange ].localePosEnd.end);
      if (contigConsistentRanges[ irange ].localePosStart.bgn -
          contigConsistentRanges[ irange ].localePosEnd.end == 0)
        fprintf( stderr, "contigConsistentRanges[ %d ].localePosStart.bgn: " F_COORD "\n",
                 irange, contigConsistentRanges[ irange ].localePosStart.bgn);
    }
  }
  fprintf( vizfile, "}\n");
}

