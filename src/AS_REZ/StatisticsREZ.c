
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
/* 	$Id: StatisticsREZ.c,v 1.10 2008-06-27 06:29:19 brianwalenz Exp $	 */

/****************************************************************************************
 *  StatisticsRez.c
 *
 *  Knut Reinert 12/99
 *
 *  contains functions to write and read statistics files for Repeat Rez
 *
 ****************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "StatisticsREZ.h"
#include "GraphCGW_T.h"
#include "ScaffoldGraph_CGW.h"
#include "DataTypesREZ.h"
#include "UtilsREZ.h"

/* ------------------------------------------------------------ */
/* functions to manipulate the GapStatisticsT struct
   This struct is the atomic piece of information that is also
   stored on disk. From it all other statistics can be computed */
/* ------------------------------------------------------------ */

#define SWITCH_THRESHOLD 2000


void init_gap_stat_struct(GapStatisticsT* g)
{
 g->lCID = g->rCID = 0;
 g->flags.all = 0;
 g->uuChunks = g->urChunks = g->ruChunks = g->rrChunks = 0;
 g->gapEstimate.mean = g->gapEstimate.variance = 0.0;
 g->gapLength.mean   = g->gapLength.variance = 0.0;
 g->exploredEdges    = g->walkedChunks = 0;
 g->uuChunks = g->urChunks = g-> ruChunks = g->rrChunks = 0;
 g->bestTooShort = g->bestTooLong = FLT_MAX;
 // ANSI C requires FLT_MAX.
}



void print_gap_stat_struct(GapStatisticsT* g, FILE* file)
{
  fprintf(file,"=-----------------------------------------------------=\n");
  fprintf(file,"=----- STATISTICS for gap contig %d <----> contig %d\n",g->lCID,g->rCID);
  fprintf(file,"=-----------------------------------------------------=\n");
  fprintf(file,"Gap walked    ? %s\n",(g->flags.bits.walked ? "YES" : "NO"));
  fprintf(file,"Gap tried     ? %s\n",(!g->flags.bits.tooShort ? "YES" : "NO"));
  fprintf(file,"Gap trivial   ? %s\n",(g->flags.bits.trivial ? "YES" : "NO"));
  fprintf(file,"Gap too short ? %s\n",(g->flags.bits.walkedTooShort ? "YES" : "NO"));
  fprintf(file,"Gap too long  ? %s\n",(g->flags.bits.walkedTooLong ? "YES" : "NO"));
  fprintf(file,"Gap maxed out ? %s\n",(g->flags.bits.maxedOut? "YES" : "NO"));

  fprintf(file,"Gap estimate (%6.2f,%6.2f)\n",g->gapEstimate.mean,sqrt(g->gapEstimate.variance));
  if( g->flags.bits.walked )
    {
      fprintf(file,"Gap length   (%6.2f,%6.2f)\n",g->gapLength.mean,sqrt(g->gapLength.variance));
      fprintf(file,"Number of edges explored = %3d / chunks walked = %3d \n",
	      g->exploredEdges,g->walkedChunks);
      fprintf(file,"Status of chunks uu = %d / ur = %d / ru = %d / rr = %d\n",g->uuChunks,g->urChunks,g->ruChunks,g->rrChunks);
    }
  else
    {
      fprintf(file,"Number of edges explored = %3d \n",
	    g->exploredEdges);

      if( g->flags.bits.walkedTooShort )
	fprintf(file,"Closest miss too short = %.2f \n",g->bestTooShort);

      if( g->flags.bits.walkedTooLong )
	fprintf(file,"Closest miss too long = %.2f \n",g->bestTooLong);
    }
}





/* ------------------------------------------------------------ */
/* functions to manipulate the ScaffoldWalkStatisticsT struct  */
/* ------------------------------------------------------------ */


/* this functions sets all builtin variables to zero.
   BEWARE that the GapStatisticsT VA is not allocated here.
   This is done either reading the gap statistics from file or
   allocating them explicitely before walking the gap */

void init_scaffold_walk_stat_struct(ScaffoldWalkStatisticsT* ws)
{
  ws->scaffoldID = NO_SCAFFOLD;
  ws->uuChunks = ws->urChunks = ws->ruChunks = ws->rrChunks = 0;
  ws->negativExploredEdges = ws->negativWalkedChunks = 0;
  ws->smallExploredEdges = ws->smallWalkedChunks = 0;
  ws->bigExploredEdges = ws->bigWalkedChunks = 0;

  ws->insertedChunks = 0;
  ws->bpsWalked = ws->bpsNotWalked = ws->bpsTried = ws->bpsMaxGap = 0;
  ws->bigGapsWalked = ws->bigGapsNotWalked = 0;
  ws->smallGapsWalked = ws->smallGapsNotWalked = ws->negativGapsWalked = 0;
  ws->negativGapsNotWalked = ws->gapsTooShort = 0;

  ws->walkedMaxedOut = ws->walkedTooShort = ws->walkedTooLong = ws->walkedTrivial = 0;
}




// Here we cycle through all gaps in a scaffold and combine their infos
void compute_scaffold_statistics(ScaffoldWalkStatisticsT* ws)
{
  int i;
  for(i=0; i< GetNumGapStatisticsTs(ws->GapStats); i++)
   {
     // get the gap stat struct from the VA
     GapStatisticsT *gapStatp = GetGapStatisticsT(ws->GapStats,i);

     if( gapStatp->flags.bits.tooShort ) // we did not try
       ws->gapsTooShort++;
     else
       {
	 if( gapStatp->flags.bits.walkedTooLong )
	   ws->walkedTooLong++;
	 if( gapStatp->flags.bits.walkedTooShort )
	   ws->walkedTooShort++;
	 if( gapStatp->flags.bits.maxedOut )
	   ws->walkedMaxedOut++;
	 if( gapStatp->flags.bits.trivial )
	   ws->walkedTrivial++;

	 // we walked it
	 if( gapStatp->flags.bits.walked  )
	   {
	     ws->uuChunks += gapStatp->uuChunks;
	     ws->urChunks += gapStatp->urChunks;
	     ws->ruChunks += gapStatp->ruChunks;
	     ws->rrChunks += gapStatp->rrChunks;
	     ws->bpsWalked += MAX(0.0,gapStatp->gapLength.mean);

	     if( gapStatp->gapLength.mean > ws->bpsMaxGap )
	       ws->bpsMaxGap = gapStatp->gapLength.mean;

	     ws->bpsWalked += gapStatp->gapLength.mean;

	     if( gapStatp->gapEstimate.mean  < 0.0 )
	       {
		 ws->negativGapsWalked++;
		 ws->negativExploredEdges +=  gapStatp->exploredEdges;
		 ws->negativWalkedChunks  +=  gapStatp->walkedChunks;

	       }
	     else
	       if( gapStatp->gapEstimate.mean < SWITCH_THRESHOLD )
		 {
		    ws->smallGapsWalked++;
		    ws->smallExploredEdges +=  gapStatp->exploredEdges;
		    ws->smallWalkedChunks  +=  gapStatp->walkedChunks;
		 }
	       else
		 {
		   ws->bigGapsWalked++;
		   ws->bigExploredEdges  +=  gapStatp->exploredEdges;
		   ws->bigWalkedChunks   +=  gapStatp->walkedChunks;
		 }
	   }
	 else // we did not walk it
	   {
	     ws->bpsNotWalked += MAX(gapStatp->gapEstimate.mean,0.0);

	     if( gapStatp->gapEstimate.mean  < 0.0 )
	       {
		 ws->negativGapsNotWalked++;
	       }
	      else
		if( gapStatp->gapLength.mean < SWITCH_THRESHOLD )
		  {
		    ws->smallGapsNotWalked++;
		  }
		else
		  {
		    ws->bigGapsNotWalked++;
		  }

	   }
       }
   }
  return;
}




void print_scaffold_walk_stat_struct(ScaffoldWalkStatisticsT* ws, FILE* file, int printGaps)
{
  if( ws->scaffoldID == NO_SCAFFOLD )
    {
      fprintf(file,"*-----Statistic struct not valid for any scaffold\n");
      return;
    }

  fprintf(file,"*-------------------------------------------------------------------------*\n");
  if( ws->scaffoldID == ALL_SCAFFOLDS )
    fprintf(file,"*----- STATISTICS for ALL scaffolds walked\n");
  else
    fprintf(file,"*----- STATISTICS for scaffold %d\n",ws->scaffoldID);

  if( printGaps )
    {
      int i;
      if( ws->scaffoldID == ALL_SCAFFOLDS )
	fprintf(file,"*----- STATISTICS for all gaps in all scaffolds\n");
      else
	fprintf(file,"*----- STATISTICS for all gaps in scaffold %d\n",ws->scaffoldID);
      for(i=0; i< GetNumGapStatisticsTs(ws->GapStats); i++)
	{
	  GapStatisticsT *gapStatp = GetGapStatisticsT(ws->GapStats,i);
	  print_gap_stat_struct(gapStatp,file);
	}

    }
  fprintf(file,"*-------------------------------------------------------------------------*\n");
  if( ws->scaffoldID == ALL_SCAFFOLDS )
    fprintf(file,"*----- Chunk statistics for all scaffolds\n");
  else
    fprintf(file,"*----- Chunk statistics for scaffold %d\n",ws->scaffoldID);
  fprintf(file,"*-------------------------------------------------------------------------*\n");
  fprintf(file,"Number of walked uu chunks        = %5d\n",ws->uuChunks);
  fprintf(file,"Number of walked ur chunks        = %5d\n",ws->urChunks);
  fprintf(file,"Number of walked ru chunks        = %5d\n",ws->ruChunks);
  fprintf(file,"Number of walked rr chunks        = %5d\n",ws->rrChunks);
  fprintf(file,"Number of inserted chunks / walked chunks / explored edges = %5d/%5d/%5d\n",
	  ws->insertedChunks,
	  ws->negativWalkedChunks+ws->smallWalkedChunks+ws->bigWalkedChunks,
	  ws->negativExploredEdges+ws->smallExploredEdges+ws->bigExploredEdges);
  fprintf(file,"*-------------------------------------------------------------------------*\n");
  if( ws->scaffoldID == ALL_SCAFFOLDS )
    fprintf(file,"*----- Chunk statistics for all scaffolds\n");
  else
    fprintf(file,"*----- Gap statistics for scaffold %d\n",ws->scaffoldID);
  fprintf(file,"*-------------------------------------------------------------------------*\n");
  {
    int negativTotal = ws->negativGapsWalked+ws->negativGapsNotWalked;
    int smallTotal   = ws->smallGapsWalked+ws->smallGapsNotWalked;
    int bigTotal     = ws->bigGapsWalked+ws->bigGapsNotWalked;
    int total        = negativTotal+smallTotal+bigTotal;

    float negativNotWalkedPercent = 0.0;
    float smallNotWalkedPercent   = 0.0;
    float bigNotWalkedPercent     = 0.0;
    float totalNotWalkedPercent   = 0.0;
    float negativWalkedPercent    = 0.0;
    float smallWalkedPercent      = 0.0;
    float bigWalkedPercent        = 0.0;
    float totalWalkedPercent      = 0.0;
    float totalWalkedBpsPercent   = 0.0;
    float totalNotWalkedBpsPercent   = 0.0;

    if( negativTotal > 0)
      negativNotWalkedPercent = 100.0 * (float) ws->negativGapsNotWalked / (float) negativTotal;

    if( smallTotal > 0)
      smallNotWalkedPercent = 100.0 *(float) ws->smallGapsNotWalked / (float) smallTotal;

    if( bigTotal > 0)
      bigNotWalkedPercent = 100.0 * (float) ws->bigGapsNotWalked / (float) bigTotal;

    if( negativTotal > 0)
      negativWalkedPercent = 100.0 * (float) ws->negativGapsWalked / (float) negativTotal;

    if( smallTotal > 0)
      smallWalkedPercent = 100.0 * (float) ws->smallGapsWalked / (float) smallTotal;

    if( bigTotal > 0)
      bigWalkedPercent = 100.0 * (float) ws->bigGapsWalked / (float) bigTotal;

    if( total > 0)
      totalWalkedPercent    = 100.0 * ( (float) ws->negativGapsWalked +
				(float) ws->smallGapsWalked +
				(float) ws->bigGapsWalked ) / (float) total;
    if( total > 0)
      totalNotWalkedPercent = 100.0 * ( (float) ws->negativGapsNotWalked +
				(float) ws->smallGapsNotWalked +
				(float) ws->bigGapsNotWalked ) / (float) total;

    if( ws->bpsNotWalked+ws->bpsWalked > 0)
      totalNotWalkedBpsPercent = 100.0 * ( (float) ws->bpsNotWalked / (float) (ws->bpsNotWalked+ws->bpsWalked));

    if( ws->bpsNotWalked+ws->bpsWalked > 0)
      totalWalkedBpsPercent = 100.0 * ( (float) ws->bpsWalked / (float) (ws->bpsNotWalked+ws->bpsWalked));

    fprintf(file,"Number of negativ     gaps walked/not walked/total = %5d/%5d/%5d\n",ws->negativGapsWalked,ws->negativGapsNotWalked,negativTotal);
    fprintf(file,"Percentage of negativ gaps walked/not walked       = %2.2f/%2.2f\n",negativWalkedPercent,negativNotWalkedPercent);
    fprintf(file,"Number of small       gaps walked/not walked/total = %5d/%5d/%5d\n",ws->smallGapsWalked,ws->smallGapsNotWalked,smallTotal);
    fprintf(file,"Percentage of small   gaps walked/not walked       = %2.2f/%2.2f\n",smallWalkedPercent,smallNotWalkedPercent);
    fprintf(file,"Number of big         gaps walked/not walked/total = %5d/%5d/%5d\n",ws->bigGapsWalked,ws->bigGapsNotWalked,bigTotal);
    fprintf(file,"Percentage of big     gaps walked/not walked       = %2.2f/%2.2f\n",bigWalkedPercent,bigNotWalkedPercent);
    fprintf(file,"Total number of       gaps walked/not walked/total = %5d/%5d/%5d\n",ws->bigGapsWalked+ws->smallGapsWalked+ws->negativGapsWalked,ws->bigGapsNotWalked+ws->smallGapsNotWalked+ws->negativGapsNotWalked,total);
    fprintf(file,"Total percentage of   gaps walked/not walked       = %2.2f/%2.2f\n",totalWalkedPercent,totalNotWalkedPercent);

    fprintf(file,"Number of gaps maxed out             = %5d/%5d\n",ws->walkedMaxedOut,total);
    fprintf(file,"Number of gaps walked but too short  = %5d/%5d\n",ws->walkedTooShort,total);
    fprintf(file,"Number of gaps walked but too long   = %5d/%5d\n",ws->walkedTooLong,total);
    fprintf(file,"Number of gaps not tried             = %5d\n",ws->gapsTooShort);
    fprintf(file,"*--------------------------------------------------------------------------*\n");
    if( ws->scaffoldID == ALL_SCAFFOLDS )
      fprintf(file,"*----- bps statistics for all scaffolds\n");
    else
      fprintf(file,"*----- bps statistics for scaffold %d\n",ws->scaffoldID);
    fprintf(file,"*--------------------------------------------------------------------------*\n");
    fprintf(file,"Number of bps walked/not walked/total   = %7d/%7d/%7d\n",
	    ws->bpsWalked,ws->bpsNotWalked,ws->bpsWalked+ws->bpsNotWalked);
    fprintf(file,"Percentage of bps walked/not walked     = %2.2f/%2.2f\n",
	    totalWalkedBpsPercent,totalNotWalkedBpsPercent);
    fprintf(file,"Number of bps tried                     = %8d\n",
	    ws->bpsTried);
    fprintf(file,"Number of bps in biggest gap            = %8d\n",
	   ws->bpsMaxGap);
    fprintf(file,"*--------------------------------------------------------------------------*\n\n");
  }
  fflush(file);
  return;
}




void allocate_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s)
{
  s->GapStats = CreateVA_GapStatisticsT(100);
  return;
}



void free_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s)
{
  DeleteVA_GapStatisticsT(s->GapStats);
  return;
}



void store_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s)
{
  char filename[200];
  FILE *output;

  AS_UTL_mkdir("stat");
  sprintf(filename,"stats/gapsInScaffold.%d.stat",s->scaffoldID);
  output = fopen(filename,"w");
  if( output == NULL )
    {
      fprintf(stderr,"=== ERROR : Could not open gap statistics file for writing in scaffold %d\n",s->scaffoldID);
      return;
    }
  else
    CopyToFileVA_GapStatisticsT(s->GapStats,output);
  fclose(output);


  sprintf(filename,"stats/scaffold.%d.stat",s->scaffoldID);
  output = fopen(filename,"w");
  if( output == NULL )
    {
      fprintf(stderr,"=== ERROR : Could not open scaffold statistics file for writing in scaffold %d\n",s->scaffoldID);
      return;
    }
  else
    {
      fprintf(output,"%d %d",s->scaffoldID,s->insertedChunks);
    }
  fclose(output);
  return;
}




int read_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s)
{

  char filename[200];
  FILE *input;

  sprintf(filename,"stats/gapsInScaffold.%d.stat",s->scaffoldID);
  input = fopen(filename,"r");
  if( input == NULL )
    {
      fprintf(stderr,"=== WARNING : Could not open gap statistics file for reading in scaffold %d\n",s->scaffoldID);
      return FALSE;
    }
  else
    s->GapStats = CreateFromFileVA_GapStatisticsT(input);
  fclose(input);

  sprintf(filename,"stats/scaffold.%d.stat",s->scaffoldID);
  input = fopen(filename,"r");
  if( input == NULL )
    {
      fprintf(stderr,"=== WARNING : Could not open scaffold statistics file for reading in scaffold %d\n",s->scaffoldID);
      return FALSE;
    }
  else
    {
      fscanf(input,"%d %d",&s->scaffoldID,&s->insertedChunks);
    }
  fclose(input);

  if( s->GapStats == NULL)
    return FALSE;
  else
    return TRUE;
}



/* ------------------------------------------------------------ */
/* functions to manipulate the WalkStatisticsT struct  */
/* ------------------------------------------------------------ */

/* the below function assumes that the accumulativ statistics for
   all scaffolds are already computed */
void compute_combined_walk_statistics(WalkStatisticsT* ws)
{
  int i;

  // now we combine all information into allScaffoldStats
  init_scaffold_walk_stat_struct(&ws->allScaffoldStats);
  ws->allScaffoldStats.scaffoldID = ALL_SCAFFOLDS;
  ws->allScaffoldStats.bpsMaxGap = 0;

  for(i=0; i<GetNumScaffoldWalkStatisticsTs(ws->ScaffoldStats); i++)
   {
     ScaffoldWalkStatisticsT *scaffStatp = GetScaffoldWalkStatisticsT(ws->ScaffoldStats,i);


     if( scaffStatp->scaffoldID != NO_SCAFFOLD)
       {
	 ws->allScaffoldStats.uuChunks += scaffStatp->uuChunks;
	 ws->allScaffoldStats.urChunks += scaffStatp->urChunks;
	 ws->allScaffoldStats.ruChunks += scaffStatp->ruChunks;
	 ws->allScaffoldStats.rrChunks += scaffStatp->rrChunks;

	 ws->allScaffoldStats.walkedTooLong  += scaffStatp->walkedTooLong;
	 ws->allScaffoldStats.walkedTooShort += scaffStatp->walkedTooShort;
	 ws->allScaffoldStats.walkedMaxedOut += scaffStatp->walkedMaxedOut;
	 ws->allScaffoldStats.walkedTrivial  += scaffStatp->walkedTrivial;

	 ws->allScaffoldStats.negativExploredEdges += scaffStatp->negativExploredEdges;
	 ws->allScaffoldStats.negativWalkedChunks  += scaffStatp->negativWalkedChunks;
	 ws->allScaffoldStats.smallExploredEdges   += scaffStatp->smallExploredEdges;
	 ws->allScaffoldStats.smallWalkedChunks    += scaffStatp->smallWalkedChunks;
	 ws->allScaffoldStats.bigExploredEdges     += scaffStatp->bigExploredEdges;
	 ws->allScaffoldStats.bigWalkedChunks      += scaffStatp->bigWalkedChunks;

	 ws->allScaffoldStats.insertedChunks += scaffStatp->insertedChunks;

	 if( scaffStatp->bpsMaxGap > ws->allScaffoldStats.bpsMaxGap )
	   ws->allScaffoldStats.bpsMaxGap = scaffStatp->bpsMaxGap;

	 ws->allScaffoldStats.bpsTried     += scaffStatp->bpsTried;
	 ws->allScaffoldStats.bpsWalked    += scaffStatp->bpsWalked;
	 ws->allScaffoldStats.bpsNotWalked += scaffStatp->bpsNotWalked;

	 ws->allScaffoldStats.bigGapsNotWalked     += scaffStatp->bigGapsNotWalked;
	 ws->allScaffoldStats.bigGapsWalked        += scaffStatp->bigGapsWalked;
	 ws->allScaffoldStats.smallGapsNotWalked   += scaffStatp->smallGapsNotWalked;
	 ws->allScaffoldStats.smallGapsWalked      += scaffStatp->smallGapsWalked;
	 ws->allScaffoldStats.negativGapsNotWalked += scaffStatp->negativGapsNotWalked;
	 ws->allScaffoldStats.negativGapsWalked    += scaffStatp->negativGapsWalked;    	        }
   }
}



/* the below function outputs a number of celagram files describing
   properties of the gaps */
void output_combined_celagram_files(WalkStatisticsT* ws)
{
  int i;
  FILE *statNumberOfNegativHops = NULL;
  FILE *statNumberOfSmallHops = NULL;
  FILE *statNumberOfBigHops = NULL;
  FILE *statNumberOfHops = NULL;
  FILE *statGapLength   = NULL;
  FILE *statGapEstimate = NULL;
  FILE *statTooShort = NULL;
  FILE *statTooLong  = NULL;

  char statFileName[256];

  // make sure that the stat directory exists
  AS_UTL_mkdir("stats");

  // now we open the celagram files
  sprintf(statFileName,"stats/number.negativ.hops.cgm");
  statNumberOfNegativHops = file_open(statFileName,"w");
  assert(NULL != statNumberOfNegativHops);
  fprintf(statNumberOfNegativHops,"Number of hops in negativ gaps\n");

  sprintf(statFileName,"stats/number.small.hops.cgm");
  statNumberOfSmallHops = file_open(statFileName,"w");
  assert(NULL != statNumberOfSmallHops);
  fprintf(statNumberOfSmallHops,"Number of hops in small gaps\n");

  sprintf(statFileName,"stats/number.big.hops.cgm");
  statNumberOfBigHops = file_open(statFileName,"w");
  assert(NULL != statNumberOfBigHops);
  fprintf(statNumberOfBigHops,"Number of hops in big gaps\n");

  sprintf(statFileName,"stats/number.hops.cgm");
  statNumberOfHops = file_open(statFileName,"w");
  assert(NULL != statNumberOfHops);
  fprintf(statNumberOfHops,"Number of hops in all gaps\n");

  sprintf(statFileName,"stats/number.gap.length.cgm");
  statGapLength = file_open(statFileName,"w");
  assert(NULL != statGapLength);
  fprintf(statGapLength,"Gap lengths in all walked gaps\n");

  sprintf(statFileName,"stats/number.gap.estimate.cgm");
  statGapEstimate = file_open(statFileName,"w");
  assert(NULL != statGapEstimate);
  fprintf(statGapEstimate,"Gap estimates in all unwalked gaps\n");

  sprintf(statFileName,"stats/number.too.short.misses.cgm");
  statTooShort = file_open(statFileName,"w");
  assert(NULL != statTooShort);
  fprintf(statTooShort,"Misses of walks were we could walk too short\n");

  sprintf(statFileName,"stats/number.too.long.misses.cgm");
  statTooLong = file_open(statFileName,"w");
  assert(NULL != statTooLong);
  fprintf(statTooLong,"Misses of walks were we could walk too long\n");

  for(i=0; i< GetNumScaffoldWalkStatisticsTs(ws->ScaffoldStats); i++)
   {
     ScaffoldWalkStatisticsT *scaffStatp = GetScaffoldWalkStatisticsT(ws->ScaffoldStats,i);

     if( scaffStatp->scaffoldID != NO_SCAFFOLD)
       {
	 int j;
	 for(j=0; j<GetNumGapStatisticsTs(scaffStatp->GapStats); j++)
	   {
	     GapStatisticsT *gapStatp = GetGapStatisticsT(scaffStatp->GapStats,j);
	     // now we write the values in the celagram files
	     // first depending on the size
	     if( gapStatp->flags.bits.walked == TRUE )
	       {
		 if( gapStatp->gapLength.mean < 0.0 )
		   {
		     fprintf(statNumberOfNegativHops,"%d ",gapStatp->walkedChunks);
		   }
		 else
		   if( gapStatp->gapLength.mean < SWITCH_THRESHOLD )
		     {
		       fprintf(statNumberOfSmallHops,"%d ",gapStatp->walkedChunks);
		     }
		   else
		     {
		       fprintf(statNumberOfBigHops,"%d ",gapStatp->walkedChunks);
		     }
		 // then general info
		 fprintf(statNumberOfHops,"%d ",gapStatp->walkedChunks);
		 fprintf(statGapLength,"%d ",(int) gapStatp->gapLength.mean);
	       }
	     else // we did not walk the guy
	       {
		 fprintf(statGapEstimate,"%d ",(int) gapStatp->gapEstimate.mean);
		 if( gapStatp->flags.bits.walkedTooShort )
		   fprintf(statTooShort,"%d ",(int) gapStatp->bestTooShort);
		 if( gapStatp->flags.bits.walkedTooLong )
		   fprintf(statTooLong,"%d ",(int) gapStatp->bestTooLong);

	       }
	   }
       }

   }

  // close all files
  fclose(statNumberOfNegativHops);
  fclose(statNumberOfSmallHops);
  fclose(statNumberOfBigHops);
  fclose(statNumberOfHops);
  fclose(statGapLength);
  fclose(statGapEstimate);
  fclose(statTooLong);
  fclose(statTooShort);
}


void allocate_walk_statistics(WalkStatisticsT *s, int num_scaff)
{
  s->ScaffoldStats = CreateVA_ScaffoldWalkStatisticsT(num_scaff);
}


/* recursively free the gap VAs in each entry of the VA ScaffoldsStats
   and then this VA itself */
void free_walk_statistics(WalkStatisticsT *ws)
{
  int i;
  for(i=0; i< GetNumScaffoldWalkStatisticsTs(ws->ScaffoldStats); i++)
    {
      ScaffoldWalkStatisticsT *scaffStatp;
      scaffStatp = GetScaffoldWalkStatisticsT(ws->ScaffoldStats,i);
      free_scaffold_walk_statistics(scaffStatp);
    }
  DeleteVA_ScaffoldWalkStatisticsT(ws->ScaffoldStats);
}

void print_all_scaffold_walk_stat_struct(char* mesg,WalkStatisticsT* ws, FILE* file, int printGaps)
{
  int y;
  for(y=0; y<GetNumScaffoldWalkStatisticsTs(ws->ScaffoldStats); y++)
    {
      ScaffoldWalkStatisticsT *scaffStatp = GetScaffoldWalkStatisticsT(ws->ScaffoldStats,y);
      fprintf(file,mesg);
      print_scaffold_walk_stat_struct(scaffStatp,file,printGaps);
    }
  return;
}
