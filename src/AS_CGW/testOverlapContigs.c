
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
static char CM_ID[] = "$Id: testOverlapContigs.c,v 1.2 2004-09-23 20:25:19 mcschatz Exp $";


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

#define AHANGSLOP 30

FILE *  File_Open (const char * Filename, const char * Mode, int exitOnFailure);
Overlap* OverlapContigs_test(NodeCGW_T *contig1, NodeCGW_T *contig2, 
                             ChunkOrientationType *overlapOrientation,
                             CDS_COORD_t minAhang, CDS_COORD_t maxAhang,
                             int computeAhang, int tryReverseComplement);
void ContigContainment_test(CIScaffoldT *scaffold,
                            NodeCGW_T *prevCI, NodeCGW_T *thisCI,
                            EdgeCGW_T *overlapEdge,
                            int tryHarder);
int straight_test(int i);

int main(int argc, char *argv[]){
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
  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("loadcgw.stderr","w");

  straight_test(1);
  // exit(1);
  
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
                                    "f:g:n:c:")) != EOF)){
      switch(ch) {
		case 'n':
		  ckptNum = atoi(argv[optind - 1]);
		  break;
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

  // cgwlog file
  sprintf(data->Output_File_Name,"%s.cgwlog_loadcgw",outputPath);
  data->logfp = File_Open (data->Output_File_Name, "w", TRUE);
  
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,
                                                  ckptNum, FALSE);    

  {
	ChunkOrientationType olapOrient;
	CDS_COORD_t minAhang, maxAhang;
        int computeAhang, tryReverseComplement;
	NodeCGW_T *contig1, *contig2;
	NodeCGW_T *scaff;

	contig1 = GetGraphNode( ScaffoldGraph->ContigGraph, 2107284);
	contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, 3460124);
	scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph,
                              contig1->scaffoldID);

	contig1->offsetAEnd.mean = 17585;
	contig1->offsetBEnd.mean = 15614;
	contig2->offsetAEnd.mean = 17473;
	contig2->offsetBEnd.mean = 15827;
	
	ContigContainment( scaff, contig1, contig2, NULL, TRUE);
  }

#if 0
  {
	ChunkOrientationType olapOrient;
	CDS_COORD_t minAhang, maxAhang;
        int computeAhang, tryReverseComplement;
	NodeCGW_T *contig1, *contig2;
	
	contig1 = GetGraphNode( ScaffoldGraph->ContigGraph, 2107284);
	contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, 3460124);
	olapOrient = AB_AB;
	computeAhang = TRUE;
	tryReverseComplement = TRUE;
	
	OverlapContigs_test( contig1, contig2, &olapOrient,
                             minAhang, maxAhang, computeAhang,
                             tryReverseComplement);
  }

  {
	ChunkOrientationType olapOrient;
	CDS_COORD_t minAhang, maxAhang;
        int computeAhang, tryReverseComplement;
	NodeCGW_T *contig1, *contig2;
	
	contig1 = GetGraphNode( ScaffoldGraph->ContigGraph, 2107284);
	contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, 3460124);
	olapOrient = BA_BA;
	computeAhang = TRUE;
	tryReverseComplement = TRUE;
	
	OverlapContigs_test( contig1, contig2, &olapOrient,
                             minAhang, maxAhang, computeAhang,
                             tryReverseComplement);
  }
#endif
}

static VA_TYPE(char) *consensus1 = NULL;
static VA_TYPE(char) *consensus2 = NULL;
static VA_TYPE(char) *quality1 = NULL;
static VA_TYPE(char) *quality2 = NULL;

Overlap* OverlapContigs_test(NodeCGW_T *contig1, NodeCGW_T *contig2, 
                             ChunkOrientationType *overlapOrientation,
                             CDS_COORD_t minAhang, CDS_COORD_t maxAhang,
                             int computeAhang,
                             int tryReverseComplement)
{
  Overlap * tempOlap1;
  int orientationFrag1, orientationFrag2;
  char *seq1, *seq2;
  double erate, thresh;
  CDS_COORD_t minlen;

  fprintf( stderr, "in OverlapContigs_test\n");

  /*
  fprintf( GlobalData->stderrc, "\ncomputing overlap for contig1: " F_CID " and contig2: " F_CID "\n", 
           contig1->id, contig2->id);
  fprintf( GlobalData->stderrc, "orientation is %c\n", (char) *overlapOrientation);
  */
  
  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;

  // if computeAhang is TRUE, allow a lot of slide in ahang
  if (computeAhang == TRUE)
  {
    minAhang = - (CDS_COORD_t) contig2->bpLength.mean;
    maxAhang = (CDS_COORD_t) contig1->bpLength.mean;
    // we subtract 3 because of an asymmetry in DP_COMPARE re AB_BA and BA_AB
    minlen -= 3;  
  }
  /* 
	 fprintf( stderr, "computeAhang is %s\n", ((computeAhang == TRUE) ? "TRUE" : "FALSE"));
	 fprintf( stderr, "minAhang is " F_COORD "\n", minAhang);
	 fprintf( stderr, "maxAhang is " F_COORD "\n", maxAhang);
  */
  if(consensus1 == NULL)
  {
    consensus1 = CreateVA_char(1024);
    consensus2 = CreateVA_char(1024);
    quality1 = CreateVA_char(1024);
    quality2 = CreateVA_char(1024);
  }else{
    ResetVA_char(consensus1);
    ResetVA_char(consensus2);
    ResetVA_char(quality1);
    ResetVA_char(quality2);
  }
  // Get the consensus sequences for both chunks from the Store
  GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
  GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);

  seq1 = Getchar(consensus1, 0);
  seq2 = Getchar(consensus2, 0);

  fprintf( stderr, "calling OverlapSequences with orient: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
		   *overlapOrientation, minAhang, maxAhang);

  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 = OverlapSequences( seq1, seq2, *overlapOrientation,
                                minAhang, maxAhang, 
                                erate, thresh, minlen, AS_FIND_ALIGN);

  if (tempOlap1 != NULL)
  {
    fprintf( stderr,
             F_CID ", " F_CID " ahang: " F_COORD ", bhang:" F_COORD "\n", 
             contig1->id, contig2->id,
             tempOlap1->begpos, tempOlap1->endpos);
    return tempOlap1;
  }
  else
  {
    fprintf( stderr, F_CID ", " F_CID " do not overlap\n",
             contig1->id, contig2->id);
    // dumpContigInfo(contig1);
    // dumpContigInfo(contig2);
    
    if ( tryReverseComplement == TRUE )
    {
      ChunkOrientationType overlapOrientationTemp = *overlapOrientation;
      
      fprintf (stderr,
               "trying reverse complementing " F_CID " and " F_CID "\n",
               contig1->id, contig2->id);
      
      Complement_Seq( seq1 );
      Complement_Seq( seq2 );
      // *overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientationTemp );
      
      fprintf( stderr, "calling OverlapSequences with orient: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
               *overlapOrientation, minAhang, maxAhang);
      
      // tempOlap1 is a static down inside of DP_Compare
      tempOlap1 = OverlapSequences( seq1, seq2, *overlapOrientation, minAhang, maxAhang, 
                                    erate, thresh, minlen, AS_FIND_ALIGN);
      
      if (tempOlap1 != NULL)
        return tempOlap1;
      else
        fprintf( stderr, F_CID ", " F_CID " do not overlap even when reverse complemented\n",
                 contig1->id, contig2->id);
    }
    else
      return NULL;	
  }
  return NULL;
}

int straight_test(int waste)
{
  char seq1[2048];  
  char seq2[2048];
  ChunkOrientationType overlapOrientation = AB_AB;
  CDS_COORD_t minAhang, maxAhang;
  double erate, thresh;
  CDS_COORD_t minlen;
  Overlap * tempOlap1;
  
  // contig 2107284
  strcpy ( seq1, "CTTCCTCTTCTTCTTCCTCTTCTTCCTTCTTCCTCTTCTCTCTTCTTCCTCTTCTTCTCCTTCTCCTCCTCCTTCTCCTCCTCCTCCTTCTTCTTCTTTTCTTCTTCTTCTTCCTATTCTTCTTCTCTTCCTCCTCCTCATTCTTCTTTTTCTGGCACAGTTGCACTACATAGCCCTGGTTAGCCTAAAACTCTCTATGTAGAACAGACTGGATTCTAAACAATTGAGATTTGCTTAATAATCTCATGAAATCCATTTGGGGATATATCTATTCATAACTATTTGGTCTATAGGACAAATTTTATAAGATAATTAATTGTATTTAGTCTGCTGGTAAATAAAAATAAATATTTAGATATAGGGGCACATTCTTGTGATCTCAGCACTTGGAAGGTAGAAGAGGGAAGATCAGGAATTCAAGGTCTAGAATTAAACTCATTGGCTATCAAGTCTAGCTAGAAGTCATGATTTTATAGTTAACACATAATTTTAAATTACTTAAAATGTAATTGAAAAACATAGATTAATGATTTGGTTTGATCTCATGCATGAATTGACAACTATGGTAAGTTTTATTGGTCTAATAGGTACATTACTATTCTTTGTAGTTTTATTTCAATTATTTTTCTTTCATATATCCTGACTACAGCTTCCCCTTCCTTCCTCTCGCCCAACCCCTTACCTCTTCTCTTTCCCAGATACACTCCTTTCCCCATGTCCTCCCATGAATATTAATTGAATACAATACAACAAGATACAATAATACCAGGCACACACCCTCACATCAAGGCTGGGCAAGACAACCAAGTAAGAGATAAAAACATCCCAAGTGCAGGGAAAAAAGACAGAAGCCCTCTTACTCCCACTATCTGAAGTCCTAAATAAATACACAAGTGTAAATATTAAACATGGATCTTTCAAATATACATTAGGAAAGTGTTCTCTTCAGGACACTTAGGAATTGCCAGGAAATGTATGGTTGTTGGGGATTGGTTCTAATGGTTTGTCTTAATTTGGCTCCCCAAATCACAGGAGAATCTGTGTGTTGAGATATTGAGAGTCCCTGTCCCCAGTTGGTTATCAATTGATTAATAAAGAGCCAACAGCCAATGACTGGGCAGGGTGAAAAATTCCCAGGAAAGAAACACAGAGGAAAGAAGGAACAGAAACACCATCATATGGAAGCAGGAAGATCAGACTTAAGAGCTATAGGATAGAAATCATCTAGCAATGAAGGTTCAAAGGGAAGCAGCCCCATGCAAGGACTGCCCAGAAGGAAACATGGCAGCAAAGATAAAATTAGATTTAGAAAGTATTAACTCAGGAATAATGAAGGGGAGTGTGTGCTAGTGGAGGGGAGATTTGGAAGTGCCCAACCACTGAGCTAGTTAAGGCAGATTAAATATAAAGGTTGCATGTGTGTCTGTTTTCATTTGAGAAAACAGTCCTTGGGTGGGTGTAGAGCCACATGCATCCACCAGGAGCATAGAGCAGCTATGAACAATTCACTGCTACATACACTAAGATCAAACTTAATAATGTAGCAAACCCTATAAAATAAGAATATATAGCAAGCCTAAAGAATAACAGAGAGAAAATGAAAACACAAATCAAGGAGTATACATGGGCTGTTTCATGACCTCAGGCACATGTATAACAGAAGACTGCCTTGTCTGCTCTCAGTAAGAGAAGCTGTACTTAATCCTGTTGAAATTTGATTCCCCAGGGTAGAGGGATGCAGGTAGGGTTGAGGTGAGGGTGGGTATGTGGGGGAGCGCTGTCTTAGAGGTGGATCAGGAATAGGGTTGATGGGGTGGAAAATTCGGGAGGAGGGACCAAGAAGGGAGTCAACTTTTGGAATGCAAATAAATAAAATAATTTAATAAAAAATACACGCATCAATTGATTGAATTCTGCTTGATCTTTTCAGTTGGTGTCACACATTGGAGGTTATGACAGGTTGGCAGAAAG");
  
// contig 3460124
  strcpy ( seq2, "CTCCTCCTCCTCCTCCTCCTCCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTTTCTTCTTCTTCTTCCTATTCTTCTTCTCTTCCTCCTCCTCATTCTTCTTTTTCTGGCACAGTTGCACTACATGGCCCTGGTTAGCCTAAAACTCTCTATGTAGACCAGACTGGATTCTAAATAATTGAGATTTGCATGGCTTTGCTTTCCAAGGTCTAGAATTAAACTCATTGGTTATCAAGTCTAGCTAGAAGTCATGATTTTATAGATAACACATAATTTTAAATTACTTAAGATGTAATTGAAAAACATAGATTAATGATTTGGTTTGATCTCATGCATGAATTGACAATTATGGTAAGTTTTATTGGTCTAATAGGTACATTATTATTCTTTGTAGTTTTATTTCAATTATTTTTCTTTCATATATCCTGACTATAGCTTCCCCTTCCTTCCTCTCCCCCAACCCCTTACCTCTTCTCTTTCCCAGATACACTCCTTTCCCCATGTCCTCCCATGAATAGTAATTGAATACAATACAACAAGATACAATAATACCAGGCACACACCCTCACATCAAGGCTGGGCAAGACAACCAAGTAAGAGATAAAAACATCCCAAGTGCAGGGAAAAGAGTCAGAGGCCCTCTTACTCCCACTATCTGGAGTCCTAAATAAATACACAAGTGTAAATATTAAACATGGATCTTTCAAATATACATTAGGAAAGTGTTCTCTTCAGGACACTTAGGAATTGCCAGGAAATGTATGGTTGTTGGGGATTGGTTCTAATGCTTTGCCTTAATTTGGCTCCCCAAATCACAGGAGAATCTGTGTGTTGAGATATTGAGAGTCCCGGTCCCCAGTTGGTTTTCAATTGATTAATAAAGAGCCAACAGCCAATGACTGGGCAGGGTGAAAAATTCCCAGGAAAGAAACAGAGGAAAGAAGGAACAGAAGCACCATCATGTGGAAGCAGGGAGATCAGACTTAAGAGCTATAGGATAGAAATCATCTAGCAATGTAGGTTCAAAGGGAAGGAGCCCCATGCGCGGACTGCCCAGAAGGAAACATGGCAGCAAAGATGAAATGAGATTTATAAAGTATTAACTCAGGAATAATGAAGGGGAGTGTGTGCTAGTGGAGGGGAGATTTGGAAGTGCCCAACCATTGAGCTAGTTAAGGCAGATTAAATATAAACGTTGCATGTATGTCTGTTTTCATTTGAGAAACCAGAGCTCTTGGGTGGGTGTAGAGCCACATGCATCCACCAGGAGCATAGAGCAGCTATGAACAATTCACTGCTACATACACTAAGATCAAACTTAATAATGTAGCAAACCCTATAAAATAAGAATATATAGCAAGCCCAAAGAATAACAGAGAGAAAATGAAAACACAAATCAAGGAGTATACATGGGCTGTTTCATGACCTCAGGCACATGTATAACAGAAGACTGCCTTGTCTGCTCTCAGTAAGAGAAGCTGTACTTAATCCTGTTGAAATTTGATTCCCCAGGGTAGAGGGATGCAGGTAGGGTTGAGGTGAGGGTGGGTATGTGGGGGAGCGCTGTCTTAGAGGTGGATCAGGAATAGGGTTGATGGGGTGGAAAACTCGGGAGGAGGGACCAAGAAGG");

  // Complement_Seq( seq0 );
  // Complement_Seq( seq1 );

  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;
  
  minAhang = - (CDS_COORD_t) strlen( seq1 );
  maxAhang = (CDS_COORD_t) strlen( seq1 );

  fprintf( stderr, "calling OverlapSequences with orient: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
           overlapOrientation, minAhang, maxAhang);

  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 = OverlapSequences( seq1, seq2, overlapOrientation,
                                minAhang, maxAhang, 
                                erate, thresh, minlen, AS_FIND_ALIGN);

  if (tempOlap1 == NULL)
	fprintf( stderr,
                 "no overlap found between seq1 and seq2 in straight_test\n");
  else
	fprintf( stderr,
                 "found overlap between seq1 and seq2 in straight_test\n");	

  return 0;
}

/****************************************************************************/
/*  
	This function is meant to be called by LeastSquares to merge
        together two contigs that have a containment relationship between them,
	but it can actually handle dovetails.  If there is no overlap,
        it asserts.
*/
void ContigContainment_test(CIScaffoldT *scaffold,
                            NodeCGW_T *prevCI, NodeCGW_T *thisCI,
                            EdgeCGW_T *overlapEdge,
                            int tryHarder)
{
  CDS_COORD_t minAhang, maxAhang;
  int flip;
  IntElementPos contigPos;
  int mergeStatus = 0;
  Overlap *contigOverlap;
  ChunkOrientationType overlapOrientation, actualOverlapOrientation;
  NodeCGW_T *leftContig, *rightContig;
  static VA_TYPE(IntElementPos) *ContigPositions = NULL;

  if(ContigPositions == NULL)
  {
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if ( min( prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean) <= 
       min( thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean))
  {
    leftContig = prevCI;
    rightContig = thisCI;
  }
  else
  {
    leftContig = thisCI;
    rightContig = prevCI;
  }

  if(ContigPositions == NULL)
  {
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if ( leftContig->offsetAEnd.mean < leftContig->offsetBEnd.mean)
  {
    // leftContig is AB
    if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) 
    {
      // rightContig is AB
      overlapOrientation = AB_AB;
      minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean -
                                leftContig->offsetAEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
    else
    {
      // rightContig is BA
      overlapOrientation = AB_BA;
      minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean -
                                leftContig->offsetAEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
  }
  else
  {
    // leftContig is BA
    if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean)
    {
      // rightContig is AB
      overlapOrientation = BA_AB;
      minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean -
                                leftContig->offsetBEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
    else 
    {
      // rightContig is BA
      overlapOrientation = BA_BA;
      minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean -
                                leftContig->offsetBEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
  }

  overlapOrientation = BA_BA;
  
  fprintf( stderr, "calling OverlapContigs with orient: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
           overlapOrientation, minAhang, maxAhang);
  contigOverlap = OverlapContigs( leftContig, rightContig,
                                  &overlapOrientation,
                                  minAhang, maxAhang, TRUE);

  if ( tryHarder )
  {
    // if no overlap found, try flipping orientation in case DPCompare asymmetry is biting us
    if (contigOverlap == NULL)
    {
      fprintf( stderr, "no overlap found between " F_CID " and " F_CID ", retrying with flipped orientation\n",
               leftContig->id, rightContig->id);
	  
      // try with the reverse orientation
      overlapOrientation =
        InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
	  
      contigOverlap = OverlapContigs( leftContig, rightContig,
                                      &overlapOrientation,
                                      minAhang, maxAhang, TRUE);
	  
      if ( contigOverlap != NULL)
      {
        CDS_COORD_t temp;
        
        temp = -contigOverlap->begpos;
        contigOverlap->begpos = -contigOverlap->endpos;
        contigOverlap->endpos = temp;
      }
      // restore the orientation
      // overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
    }
	
    // if still no overlap found, try maxing out hangs
    if (contigOverlap == NULL)
    {
      CDS_COORD_t maxLength;
      
      fprintf( stderr, "no overlap found between " F_CID " and " F_CID ", retrying with max AHANGSLOP\n",
               leftContig->id, rightContig->id);
	  
      maxLength = max( leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
      fprintf( stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
               (char) overlapOrientation, -maxLength, maxLength);
	  
      contigOverlap = OverlapContigs( leftContig, rightContig,
                                      &overlapOrientation,
                                      -maxLength, maxLength, FALSE);	  
    }
  
    // if still no overlap found, try flipping orientation and maxing out hangs
    if (contigOverlap == NULL)
    {
      CDS_COORD_t maxLength;
      
      fprintf( stderr, 
               "no overlap found between " F_CID " and " F_CID ", retrying with flipped orientation and max AHANGSLOP\n",
               leftContig->id, rightContig->id);
	  
      maxLength = max( leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
      // try with the reverse orientation
      overlapOrientation =
        InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
	  
      fprintf( stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
               (char) overlapOrientation, -maxLength, maxLength);
	  
      contigOverlap = OverlapContigs( leftContig, rightContig,
                                      &overlapOrientation,
                                      -maxLength, maxLength, FALSE);
	  
      if ( contigOverlap != NULL)
      {
        CDS_COORD_t temp;
        
        temp = -contigOverlap->begpos;
        contigOverlap->begpos = -contigOverlap->endpos;
        contigOverlap->endpos = temp;
      }
	  
      // restore the orientation
      // overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
    }
  }

  if (contigOverlap == NULL)
  {
    fprintf( stderr, "no overlap found between " F_CID " and " F_CID ", aborting...\n",
             leftContig->id, rightContig->id);
    dumpContigInfo(leftContig);
    dumpContigInfo(rightContig);
    assert(0);
  }	
  
  if (contigOverlap->begpos < 0) // contigs need to be reversed
  {
    leftContig = thisCI;
    rightContig = prevCI;
    // adjust Overlap fields for later use in positioning
    contigOverlap->begpos = - contigOverlap->begpos;
    contigOverlap->endpos = - contigOverlap->endpos;	
    
    switch(overlapOrientation)
    {
      case AB_AB:
      case BA_BA:
        actualOverlapOrientation = overlapOrientation;	  
        break;
      case AB_BA:
        actualOverlapOrientation = BA_AB;
        break;
      case BA_AB:
        actualOverlapOrientation = AB_BA;
        break;
      default:
        assert(0);
    }
    fprintf( stderr, "* Switched right-left  orientation went from %c to %c\n",
             overlapOrientation, actualOverlapOrientation);
  }
  else
  {
    actualOverlapOrientation = overlapOrientation;
  }
  
  fprintf( stderr, "* Containing contig is " F_CID " contained contig is " F_CID " ahg:" F_COORD "  bhg:" F_COORD " orient:%c\n",
           leftContig->id, rightContig->id,
           contigOverlap->begpos, contigOverlap->endpos,
           actualOverlapOrientation);

  fprintf( stderr,
           "* Initial Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
           leftContig->id, 
           leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
           rightContig->id, 
           rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
  
  // assume we leave the leftContig where it is
  if ( actualOverlapOrientation == AB_AB)
  {
    rightContig->offsetAEnd.mean =
      leftContig->offsetAEnd.mean + contigOverlap->begpos;
    rightContig->offsetBEnd.mean =
      rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == AB_BA)
  {
    rightContig->offsetBEnd.mean =
      leftContig->offsetAEnd.mean + contigOverlap->begpos;
    rightContig->offsetAEnd.mean =
      rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == BA_AB)
  {
    rightContig->offsetAEnd.mean =
      leftContig->offsetBEnd.mean + contigOverlap->begpos;
    rightContig->offsetBEnd.mean =
      rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  if ( actualOverlapOrientation == BA_BA)
  {
    rightContig->offsetBEnd.mean =
      leftContig->offsetBEnd.mean + contigOverlap->begpos;
    rightContig->offsetAEnd.mean =
      rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }
  
  contigPos.ident = leftContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = leftContig->offsetAEnd.mean;
  contigPos.position.end = leftContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  contigPos.ident = rightContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = rightContig->offsetAEnd.mean;
  contigPos.position.end = rightContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  
  fprintf( stderr,
           "* Final Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
           leftContig->id, 
           leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
           rightContig->id, 
           rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
  
  flip = (leftContig->offsetBEnd.mean < leftContig->offsetAEnd.mean);
  if(flip)
  {
    mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions,
                                          leftContig->offsetBEnd,
                                          leftContig->offsetAEnd);
  }
  else
  {
    mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions,
                                          leftContig->offsetAEnd,
                                          leftContig->offsetBEnd);
  }
  assert(mergeStatus == TRUE);
}

