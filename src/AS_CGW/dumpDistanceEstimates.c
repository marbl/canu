
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

static char CM_ID[] = "$Id: dumpDistanceEstimates.c,v 1.20 2007-03-27 07:31:59 brianwalenz Exp $";

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
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"
#include "AS_ALN_aligners.h"
#include "AS_ALN_forcns.h"

#include <SYS_UIDcommon.h>
#include <SYS_UIDclient.h>

//  Writes, to stderr, information about the library distances.  Also
//  creates a stat/ directory with stats.

//  params for which operations to do

//  param for min num samples to use

//  how to detect bad distributions?



//  This is appropriate for getting exactly one UID.  If you want more
//  than one, see SYS_UIDclient_simple.c
//
CDS_UID_t
dde_getUID(int real) {
  CDS_UID_t  uid         = 0;
  CDS_UID_t  interval[4] = {0};
  int        sta         = 0;

  get_uids(1, interval, real);
  sta = get_next_uid(&uid, real);
  if (sta != UID_CODE_OK) {
    fprintf(stderr, "Couldn't get a UID (get_uids()).\n");
    exit(1);
  }

  return(uid);
}




int
dde_compareDistances( const void *s1, const void *s2) {
  const MateInfoT * t1 = (const MateInfoT *) s1;
  const MateInfoT * t2 = (const MateInfoT *) s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if (t1->samples < t2->samples)
    return -1;
  else if (t1->samples > t2->samples)
    return 1;
  else 
    return 0;
}


int
dde_compareIntegers(const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}


void
dde_stats(int         operateOnNodes,
          int32       minSamplesForOverride,
          char       *instance_label,
          CDS_UID_t  *libMap) {

  GraphCGW_T *graph = NULL;
  GraphNodeIterator nodes;
  NodeCGW_T *node;
  DistT *dptr;
  int i, ii, maxSamples = 1;
  MultiAlignT *ma = CreateEmptyMultiAlignT();

  int numPotentialRocks = 0;
  int numPotentialStones = 0;

  int NN = GetNumDistTs(ScaffoldGraph->Dists);

  VA_TYPE(CDS_CID_t) *dptrFrags[NN];
  VA_TYPE(CDS_CID_t) *dptrMates[NN];
  

  if (operateOnNodes == UNITIG_OPERATIONS)
    graph = ScaffoldGraph->CIGraph;
  else if (operateOnNodes == CONTIG_OPERATIONS)
    graph = ScaffoldGraph->ContigGraph;
  else if (operateOnNodes == SCAFFOLD_OPERATIONS)
    graph = ScaffoldGraph->CIGraph;
  
  // Initialize computed mean/variance and counters, and bounds
  for(i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dptr = GetDistT(ScaffoldGraph->Dists,i);

    dptr->mu             = 0.0;
    dptr->sigma          = 0.0;
    dptr->numSamples     = 0.0;
    dptr->numReferences  = 0;
    dptr->min            = CDS_COORD_MAX;
    dptr->max            = CDS_COORD_MIN;
    dptr->bsize          = 0;
    dptr->numBad         = 0;
    dptr->histogram      = NULL;

    // Lower and upper are based on nominal mean and stddev
    dptr->lower    = dptr->mean - CGW_CUTOFF*dptr->stddev;
    dptr->upper    = dptr->mean + CGW_CUTOFF*dptr->stddev;
      
    // starting from scratch or start from a checkpoint
    //
    if(dptr->samples)
      ResetVA_CDS_COORD_t(dptr->samples);
    else
      dptr->samples = CreateVA_CDS_COORD_t(16);
    
    dptrFrags[i] = CreateVA_CDS_CID_t(256);
    dptrMates[i] = CreateVA_CDS_CID_t(256);
  }
  

  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);

  
  while(NULL != (node = NextGraphNodeIterator(&nodes))) {
    CDS_CID_t i;
    int numFrags;
    int numExternalLinks = 0;
    
    // Don't waste time loading singletons for this
    if(node->flags.bits.isChaff)
      continue;
    
    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, node->id, graph->type == CI_GRAPH);
    numFrags  = GetNumIntMultiPoss(ma->f_list);
    
    if (numFrags < 2)
      continue;
    
    for( i = 0; i < numFrags; i++) {
      IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
      CIFragT *frag, *mate;
      CDS_COORD_t dist;
      
      frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mp->sourceInt);
      assert(frag->iid == mp->ident);

      // This is important for keeping our computation as local as possible.
      // We skip fragments that have external links only, or no links
      //
      if (frag->flags.bits.hasMate == 0)
        continue;

      // the typical case
      if ((operateOnNodes == CONTIG_OPERATIONS) &&
          (frag->flags.bits.hasInternalOnlyContigLinks == 0))
        continue;
      
      mate = GetCIFragT(ScaffoldGraph->CIFrags,frag->mateOf);
      
      if ((operateOnNodes == UNITIG_OPERATIONS) && (mate != NULL) && (mate->cid != frag->cid))
        numExternalLinks++;

      
      // If the id of the current fragment is greater than its mate, skip it, to
      // avoid double counting.
      //
      if ((CDS_CID_t)mp->sourceInt > frag->mateOf)
        continue;
      

      // now make sure the 5p end is less than the 3p end
      // what about outtie mates???
      //
      if ( frag->offset5p.mean > frag->offset3p.mean)
        continue;
      
      assert(mate != NULL && mate->mateOf == (CDS_CID_t)mp->sourceInt);

      dptr = GetDistT(ScaffoldGraph->Dists, frag->dist);
      dptr->numReferences++;
      

      if (operateOnNodes == UNITIG_OPERATIONS) {
        NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);

        if (frag->cid != mate->cid)
          continue;
        
        dist = mate->offset5p.mean - frag->offset5p.mean;
        
        if(	getCIFragOrient(mate) == getCIFragOrient(frag)) {
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dptr->numBad++;
          continue;
        }

        // try to sample fairly by only doing mates where ones of any length could live
        if ( frag->offset5p.mean + dptr->mean + 5 * dptr->stddev > unitig->bpLength.mean)
          continue;

        if((frag->flags.bits.innieMate && getCIFragOrient(frag) == B_A) ||
           (!frag->flags.bits.innieMate && getCIFragOrient(frag) == A_B) )
          dist = -dist;
      } else if (operateOnNodes == CONTIG_OPERATIONS) {
        ContigT *contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
        
        assert(frag->contigID == mate->contigID);
        if(GetContigFragOrient(mate) == GetContigFragOrient(frag)) {
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dptr->numBad++;
          continue;
        }
        
        // try to sample fairly by only doing mates where ones of any length could live
        if ( frag->contigOffset5p.mean + dptr->mean + 5 * dptr->stddev > contig->bpLength.mean)
          continue;
        
        dist =  mate->contigOffset5p.mean - frag->contigOffset5p.mean; 
        //   -------------------->          <----------------------
        //     mate                 innie            frag 
        
        //   <-------------------           ---------------------->
        //     mate                 outie            frag 
        
        if((frag->flags.bits.innieMate && GetContigFragOrient(frag) == B_A) ||
           (!frag->flags.bits.innieMate && GetContigFragOrient(frag) == A_B) )
          dist =  -dist;
        
      } else if (operateOnNodes == SCAFFOLD_OPERATIONS) {
        NodeCGW_T *fragContig, *mateContig;
        CDS_COORD_t fragLeftEnd, fragRightEnd;
        CDS_COORD_t mateLeftEnd, mateRightEnd;
        int mateScaffoldOrientation, fragScaffoldOrientation;
        
        fragContig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
        AssertPtr(fragContig);
        
        mateContig = GetGraphNode( ScaffoldGraph->ContigGraph, mate->contigID);
        AssertPtr(mateContig);
        
        
        // we want them to be in the same scaffold
        if ( fragContig->scaffoldID != mateContig->scaffoldID || fragContig->scaffoldID == -1)
          continue;
        
        GetFragmentPositionInScaffold( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation);
        GetFragmentPositionInScaffold( mate, &mateLeftEnd, &mateRightEnd, &mateScaffoldOrientation);
        
        // try to sample fairly by only doing mates where ones of any length could live
        {
          NodeCGW_T *scaff, *extremeContig;
          CDS_COORD_t contigLeftEnd, contigRightEnd;
          int contigScaffoldOrientation;
          
          // grab scaffold
          scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, fragContig->scaffoldID);
          extremeContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.BEndCI);
          GetContigPositionInScaffold ( extremeContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);
          
          if ( fragLeftEnd + dptr->mean + 5 * dptr->stddev > contigRightEnd)
            continue;
        }
        
        
        if (fragScaffoldOrientation == mateScaffoldOrientation) {
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dptr->numBad++;
          continue;
        }
        
        if (frag->flags.bits.innieMate) {
          if (fragScaffoldOrientation == 0) // frag ---->  <---- mate
            dist = mateRightEnd - fragLeftEnd;		  
          else                              // mate ---->  <---- frag
            dist = fragRightEnd - mateLeftEnd;
        } else {  //  outtie pair
          if (fragScaffoldOrientation == 0) // mate <----  ----> frag
            dist = fragRightEnd - mateLeftEnd;		  
          else                              // frag <----  ----> mate
            dist = mateRightEnd - fragLeftEnd;
        }
        
        //if (dist < 0)
        //  fprintf( stderr, "frag, mate: " F_CID ", " F_CID " have negative dist: " F_COORD "\n",
        //           frag->iid, mate->iid, dist);
      }





      
      if (dist > 0 && dist < dptr->min)
        dptr->min = dist;
      if (dist > dptr->max)
        dptr->max = dist;
      
      AppendCDS_COORD_t( dptr->samples, &dist);
      AppendCDS_CID_t( dptrFrags[frag->dist], &frag->iid);
      AppendCDS_CID_t( dptrMates[frag->dist], &mate->iid);
      
     
      // See if the mate distance implied is outside of a 5-sigma range
        
      if (dist < dptr->lower || dist > dptr->upper) {
        //if (dist > 5000)
        // fprintf(stderr,"* (" F_CID "," F_CID ") lib:" F_CID " is bad due to distance problems " F_COORD " is outside [" F_COORD "," F_COORD "]\n",
        //	frag->iid, mate->iid, frag->dist, dist, dptr->lower, dptr->upper);
        frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
        mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
        dptr->numBad++;


      } else {
        frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
        mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          
        //if (dist > 5000)
        //fprintf(stderr,"* (" F_CID "," F_CID ") lib:" F_CID " is trusted due to distance " F_COORD " being inside [" F_COORD "," F_COORD "]\n",
        //	frag->iid, mate->iid, frag->dist, dist, dptr->lower, dptr->upper);
          
        dptr->numSamples++;
        dptr->mu += dist;
        dptr->sigma += (((double)dist)*(double)dist);

      }
    }

    // Mark unitigs as potential Rocks and Stones
    if (operateOnNodes == UNITIG_OPERATIONS) {
      int rock = FALSE;
      int stone = FALSE;
      switch(numExternalLinks){
        case 0:
          stone = TRUE;
          rock = FALSE;
          break;
        case 1:
          stone = TRUE;
          rock = FALSE;
          numPotentialStones++;
          break;
        case 2:
        default:
          stone = rock = TRUE;
          numPotentialRocks++;
          break;
      }
      node->flags.bits.isPotentialRock = rock;
      node->flags.bits.isPotentialStone = stone;
    }
  }
  
  if (operateOnNodes == UNITIG_OPERATIONS) {
    fprintf(stderr,
            "* ComputeMatePairStats has marked %d/%d unitigs as potential rocks +  %d/%d as potential stones\n",
            numPotentialRocks, (int) GetNumGraphNodes(graph),
            numPotentialStones, (int) GetNumGraphNodes(graph));
  }
  
  
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    if ( GetNumCDS_COORD_ts( GetDistT(ScaffoldGraph->Dists, i)->samples) > maxSamples)
      maxSamples = GetNumCDS_COORD_ts( GetDistT(ScaffoldGraph->Dists, i)->samples);
  }
  
  fprintf( stderr, "maxSamples = %d\n", maxSamples);
  
  // now sort the samples, mates, and frags arrays, based on samples
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    // MateInfoT *matePairs;
    MateInfoT matePairs[maxSamples]; // matePairs[maxSamples];
    int icnt;
    CDS_COORD_t newLower, newUpper;
    CDS_COORD_t median = 0, lowerSigma = 0, upperSigma = 0;
    
    dptr = GetDistT(ScaffoldGraph->Dists, i);
    if(dptr->numReferences == 0)
      continue;
    if (dptr->numSamples == 0 || dptr->numSamples == 1)
      continue;
    
    for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++) {
      matePairs[icnt].samples = *GetCDS_COORD_t( dptr->samples, icnt);
      matePairs[icnt].frags = *GetCDS_CID_t( dptrFrags[i], icnt);
      matePairs[icnt].mates = *GetCDS_CID_t( dptrMates[i], icnt);
    }
    
    qsort( matePairs, GetNumCDS_COORD_ts( dptr->samples ),
           sizeof(MateInfoT), &dde_compareDistances);

    // now find median
    median = matePairs[ (int) (GetNumCDS_COORD_ts( dptr->samples ) / 2)].samples;
    
    // find lower and upper "sigmas"
    lowerSigma =   median - matePairs[ (int) ((0.5 - 0.34) * GetNumCDS_COORD_ts( dptr->samples ))].samples;
    upperSigma = - median + matePairs[ (int) ((0.5 + 0.34) * GetNumCDS_COORD_ts( dptr->samples ))].samples;	

    newLower = median - 5 * MAX (lowerSigma, upperSigma);
    if ( newLower < 0 )
      newLower = 0;
    newUpper = median + 5 * MAX (lowerSigma, upperSigma);
    
#if 0
    fprintf( stderr, "\nlib " F_CID ", numSamples: %d, orig mean, sig: ( %.2f, %.2f), calc mean, sig: (%.2f, %.2f) median: " F_COORD "\n",
             i, dptr->numSamples, dptr->mean, dptr->stddev, dptr->mu/dptr->numSamples, 
             sqrt((dptr->sigma - (dptr->mu*dptr->mu)/dptr->numSamples) / (dptr->numSamples - 1)),
             median);
    fprintf( stderr, "dptr->lower: " F_COORD ", dptr->upper: " F_COORD "\n", dptr->lower, dptr->upper);
    fprintf( stderr, "  dptr->min: " F_COORD ",   dptr->max: " F_COORD "\n", dptr->min, dptr->max);
    fprintf( stderr, "lowerSigma: " F_COORD ", upperSigma: " F_COORD "\n", lowerSigma, upperSigma);						  
    fprintf( stderr, "newLower: " F_COORD ", newUpper: " F_COORD "\n", newLower, newUpper);
#endif
    
    // now reset the trusted flag if necessary
    // first see if there are edges marked untrusted that are now considered trusted
    // lower set
    if (dptr->lower > newLower) {
      for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++) {
        if (matePairs[icnt].samples < dptr->lower && matePairs[icnt].samples > newLower) {
          CIFragT *frag, *mate;
          
          frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
          mate = GetCIFragT( ScaffoldGraph->CIFrags,
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
          //fprintf( stderr, "1 reclassifying samples[%d] (" F_CID ", " F_CID ") from UNTRUSTED to TRUSTED\n", 
          //         icnt, frag->iid, mate->iid);

          frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          dptr->numBad--;
          dptr->numSamples++;
          dptr->mu += matePairs[icnt].samples;
          dptr->sigma += (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
        }
      }
    }
    
    // upper set
    if (dptr->upper < newUpper) {
      for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++) {	
        if (matePairs[icnt].samples > dptr->upper && matePairs[icnt].samples < newUpper) {
          CIFragT *frag, *mate;
          
          frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
          mate = GetCIFragT( ScaffoldGraph->CIFrags,
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
          //fprintf( stderr, "2 reclassifying samples[%d] (" F_CID ", " F_CID ") from UNTRUSTED to TRUSTED\n", 
          //         icnt, frag->iid, mate->iid);
          
          frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          dptr->numBad--;
          dptr->numSamples++;
          dptr->mu += matePairs[icnt].samples;
          dptr->sigma += (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
        }
      }
    }



    
    // now see if there are edges marked trusted that are now considered untrusted
    // lower set
    if (dptr->lower < newLower) {
      for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++) {
        if (matePairs[icnt].samples > dptr->lower && matePairs[icnt].samples < newLower) {
          CIFragT *frag, *mate;
          
          frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
          mate = GetCIFragT( ScaffoldGraph->CIFrags,
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
          //fprintf( stderr, "3 reclassifying samples[%d] (" F_CID ", " F_CID ") from TRUSTED to UNTRUSTED\n", 
          //         icnt, frag->iid, mate->iid);
          
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dptr->numBad++;
          dptr->numSamples--;
          dptr->mu -= matePairs[icnt].samples;
          dptr->sigma -= (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
        }
      }
    }
    // upper set
    if (dptr->upper > newUpper) {
      for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++) {	
        if (matePairs[icnt].samples < dptr->upper && matePairs[icnt].samples > newUpper) {
          CIFragT *frag, *mate;
          
          frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
          mate = GetCIFragT( ScaffoldGraph->CIFrags,
                             GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
          //fprintf( stderr, "4 reclassifying samples[%d] (" F_CID ", " F_CID ") from TRUSTED to UNTRUSTED\n", 
          //         icnt, frag->iid, mate->iid);
          
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dptr->numBad++;
          dptr->numSamples--;
          dptr->mu -= matePairs[icnt].samples;
          dptr->sigma -= (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
        }
      }
    }
  }



  
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DeleteVA_CDS_CID_t( dptrFrags[i] );
    DeleteVA_CDS_CID_t( dptrMates[i] );
  }
  


  
  /* now set mean, stddev, size & number of buckets */
  for(i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    dptr = GetDistT(ScaffoldGraph->Dists, i);
    if(dptr->numReferences == 0) {
      //fprintf(stderr,"* No references for mate distance record %d\n", i);
      continue;
    }
    if (dptr->numSamples == 0) {
      //fprintf(stderr,"* No samples for mate distance record %d... numBad = %d/%d\n", 
      //        i, dptr->numBad, dptr->numReferences);
      continue;
    }
    if(dptr->numSamples == 1)
      dptr->sigma = -1.0;
    else
      dptr->sigma = sqrt((dptr->sigma - (dptr->mu*dptr->mu)/dptr->numSamples) /
                         (dptr->numSamples - 1));
    
    dptr->mu = dptr->mu/dptr->numSamples;
    
    if(dptr->numSamples > minSamplesForOverride) {
      fprintf(stderr, "Distance %3d mean %9.2f (nominal %9.2f) stddev %9.2f (nominal %9.2f) numSamples %6d numRef %6d\n",
              i, 
              dptr->mu, dptr->mean,
              dptr->sigma, dptr->stddev,
              dptr->numSamples, dptr->numReferences);

      if (dptr->mu - 3.0 * dptr->sigma < 0.0) {
        fprintf(stderr, "Distance %3d is INVALID!  Possibly bimodal?\n", i);
      } else {
        fprintf(stdout, "{DST\n");
        fprintf(stdout, "act:R\n");
        fprintf(stdout, "acc:"F_UID"\n", libMap[i]);
        fprintf(stdout, "mea:%f\n", dptr->mu);
        fprintf(stdout, "std:%f\n", dptr->sigma);
        fprintf(stdout, "}\n");
      }

      //  We don't update here, just report
      //dptr->mean   = dptr->mu;
      //dptr->stddev = dptr->sigma;
    }


    dptr->bsize = dptr->sigma/CGW_NUM_BUCKETS;
    if (dptr->bsize > 1.) {
      dptr->bnum = (dptr->max-dptr->min)/dptr->bsize + 1;
      dptr->bsize = ((float)(dptr->max-dptr->min))/dptr->bnum;
    } else {
      dptr->bnum = 1;
      dptr->bsize = dptr->max-dptr->min;
    }

    safe_free(dptr->histogram);
    dptr->histogram = (int32 *) safe_malloc(sizeof(int32)*dptr->bnum);

    {
      int	j;
      for (j=0; j < dptr->bnum; ++j)
        dptr->histogram[j] = 0;
    }
  }

  
  // Now bucketize data
  {
    int32 numDists = GetNumDistTs(ScaffoldGraph->Dists);
    DistT *dist = GetDistT(ScaffoldGraph->Dists,1);

    for(ii = 1; ii < numDists ; ii++, dist++){
      int numSamples = GetNumCDS_COORD_ts(dist->samples);
      CDS_COORD_t *samplep = GetCDS_COORD_t(dist->samples,0);
      int j;

      if(!dist->samples || dist->bsize == 0.0)
        continue;

      for(j = 0; j < numSamples ; j++, samplep++) {
        CDS_COORD_t sample = *samplep;
        int32 binNum = ((float)sample - (float)dist->min)/dist->bsize;
        binNum = MIN(binNum, dist->bnum -1);
        binNum = MAX(binNum,0); // sample can be negative
        dist->histogram[binNum]++;
      }
    }
  }



  //  output a histogram file for each library
  //
  for( ii = 1; ii < GetNumDistTs(ScaffoldGraph->Dists); ii++) {
    int j;
    DistT *dist = GetDistT(ScaffoldGraph->Dists,ii);
    if(dist->samples && dist->bsize > 0.0 && dist->numSamples > 0) {
      FILE        *fout;
      char         filename[1024];
      int          numSamples          = GetNumCDS_COORD_ts(dist->samples);
      CDS_COORD_t  samps[numSamples];
      CDS_COORD_t *samplep             = 0L;

#if 0
      fprintf(stderr,"* Distance Record %d min:" F_COORD " max:" F_COORD "  mu:%g sigma:%g samples:%d bad:%d references:%d\n",
              ii, dist->min, dist->max, dist->mu, dist->sigma, dist->numSamples, dist->numBad, dist->numReferences);
      for(j = 0; j < dist->bnum; j++)
        fprintf(stderr,"* [%5d,%5d]\t%d\n",
                (int32)(dist->min + j * dist->bsize),
                (int32)(dist->min + (j + 1) * dist->bsize),
                dist->histogram[j]);
#endif

      sprintf( filename, "%s.distlib_%03d.cgm", instance_label, ii);
      fout = fopen( filename,"w");
      if (fout == NULL)
        fprintf(stderr, "Failed to open '%s' for writing.\n", filename), exit(1);

      fprintf(fout, "lib %d mu %g sigma %g\n", ii, dist->mu, dist->sigma );

      samplep = GetCDS_COORD_t(dist->samples,0);
      for(j = 0; j < numSamples; j++, samplep++)
        samps[j] = *samplep;

      qsort(samps, numSamples, sizeof(CDS_COORD_t), &dde_compareIntegers);

      for (j = 0; j < numSamples; j++)
        fprintf(fout, "%d\n", samps[j]);

      fclose(fout);
    }
  }
  
  DeleteMultiAlignT(ma);
}





//  Swiped from dumpGatekeeper
//
CDS_UID_t*
buildLibraryIIDmap(char *gkpStoreName) {
  GateKeeperStore         *gkpStore;
  GateKeeperLibraryRecord  gkpl;
  StoreStat                stat;
  int64                    i;
  CDS_UID_t               *map;

  gkpStore = openGateKeeperStore(gkpStoreName, FALSE);

  statsStore(gkpStore->lib, &stat);
  fprintf(stderr,"* Stats for Dist Store are first:" F_S64 " last :" F_S64 "\n",
          stat.firstElem, stat.lastElem);

  map = (CDS_UID_t *)safe_malloc(sizeof(CDS_UID_t) * (stat.lastElem+1));

  for (i=0; i<=stat.lastElem; i++)
    map[i] = 0;

  for (i=1; i<=stat.lastElem; i++) {
    getGateKeeperLibraryStore(gkpStore->lib, i, &gkpl);
    map[i] = gkpl.libraryUID;

    fprintf(stderr,"* Dist " F_S64 " UID:" F_UID " del:%d red:%d mean:%f std:%f batch(" F_U16 "," F_U16 ") prevID:" F_IID " prevInstanceID: " F_IID "\n",
            i,
            gkpl.libraryUID,
            gkpl.deleted,
            gkpl.redefined,
            gkpl.mean,
            gkpl.stddev,
            gkpl.birthBatch,
            gkpl.deathBatch,
            gkpl.prevID,
            gkpl.prevInstanceID);
  }

  closeGateKeeperStore(gkpStore);

  return(map);
}




int
main( int argc, char **argv) {
  Global_CGW  *data                  = 0L;
  int          ckptNum               = NULLINDEX;
  int          minSamplesForOverride = 50;  //  Was 1000, but that doesn't work for smaller assemblies.
  CDS_UID_t   *libMap                = 0L;
  int          realUIDs              = FALSE;
  CDS_UID_t    uidStart              = 987654321987654321ULL;
  
  int arg=1;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc   = stderr;
  data->timefp    = stderr;

  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      arg++;
      strcpy(data->Gatekeeper_Store_Name, argv[arg]);
    } else if (strncmp(argv[arg], "-prefix", 2) == 0) {
      arg++;
      strcpy(data->File_Name_Prefix, argv[arg]);
    } else if (strncmp(argv[arg], "-nckp", 2) == 0) {
      arg++;
      ckptNum = atoi(argv[arg]);

    } else if (strncmp(argv[arg], "-u", 2) == 0) {
      realUIDs = TRUE;
    } else if (strncmp(argv[arg], "-s", 2) == 0) {
      realUIDs = FALSE;
      set_start_uid(STR_TO_UID(argv[++arg], NULL, 10));
#ifdef USE_SOAP_UID
    } else if (strncmp(argv[arg], "-E", 2) == 0) {
      SYS_UIDset_euid_server(argv[++arg]);
#endif
    } else {
      fprintf(stderr, "%s: unknown arg %s\n", argv[0], argv[arg]);
    }

    arg++;
  }

  if (ckptNum == NULLINDEX) {
    fprintf(stderr, "usage: %s [opts] -gkp gkpStore -prefix asmprefix -nckp ckptNum\n", argv[0]);
    fprintf(stderr, "  -u                 Use real UIDs from the UID server\n");
    fprintf(stderr, "  -s uid             Use the fake UID 'uid', default uid=987654321987654321\n");
    fprintf(stderr, "  -E <server:port>   Use EUID server server:port insteat of tools.tigr.org:8190 (does NOT imply -u)\n");
    exit(1);
  }

  //  We need to map the library IID to UID, so we can create a gatekeeper update.
  //  Alternatively, we could modify the gatekeeperstore directly here.
  //
  libMap = buildLibraryIIDmap(data->Gatekeeper_Store_Name);

  //  LoadScaffoldGraphFromCheckpoint wants to CheckCIScaffoldT()
  //  which can RecomputeOffsetsInScaffold(), which can eventually,
  //  try to get an overlap.  Unless this is set, it bombs.
  //
  GlobalData->aligner=Local_Overlap_AS_forCNS;

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix, ckptNum, FALSE);

  fprintf(stdout, "{BAT\n");
  fprintf(stdout, "bna:Distance Update\n");
  fprintf(stdout, "crt:0\n");
  fprintf(stdout, "acc:"F_UID"\n", dde_getUID(realUIDs));
  fprintf(stdout, "com:\n");
  fprintf(stdout, "Generated by %s\n", argv[0]);
  fprintf(stdout, ".\n");
  fprintf(stdout, "}\n");

  fprintf(stderr, "================================================================================\n");
  fprintf(stderr, "= SCAFFOLDS                                                                    =\n");
  fprintf(stderr, "================================================================================\n");
  dde_stats(SCAFFOLD_OPERATIONS, minSamplesForOverride, "scaffold_final", libMap);

  fprintf(stderr, "================================================================================\n");
  fprintf(stderr, "= CONTIGS                                                                      =\n");
  fprintf(stderr, "================================================================================\n");
  dde_stats(CONTIG_OPERATIONS, minSamplesForOverride, "contig_final", libMap);

  //  eCR should call CheckCIScaffoldTs() before dumping a checkpoint -- so we won't have to do any
  //  thing here when we load the checkpoint (change the TRUE to a FALSE then).
  //
  //CheckpointScaffoldGraph(ScaffoldGraph, 1);

  DestroyScaffoldGraph(ScaffoldGraph);
  DeleteGlobal_CGW(data);

  exit(0);
}
