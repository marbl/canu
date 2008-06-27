
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
/* 	$Id: InputDataTypes_CGW.h,v 1.16 2008-06-27 06:29:14 brianwalenz Exp $	 */
/****************************************************************************
 *  InputDataTypes_CGW
 *
 *  Saul A. Kravitz 9/99
 *
 *  Definitions for the data structures used to hold fragments and distance records
 *
 ***************************************************************************/
#ifndef INPUTDATATYPES_CGW_H
#define INPUTDATATYPES_CGW_H

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

typedef struct {
  double  mean;
  double variance;
}LengthT;

// Do the arithmetic and stats on two unordered pairs of LengthTs
// If resulting variance is negative assert
void ComputeIntervalLength(LengthT *result,
			   LengthT *aEndA, LengthT *bEndA,
			   LengthT *aEndB, LengthT *bEndB);


// Do the arithmetic and stats on a pair of LengthTs
// If resulting variance is negative assert
void ComputeLength(LengthT *result,
		   LengthT *length1, LengthT *length2);


typedef enum {X_X = 'X', A_B = 'F', B_A = 'R'} FragOrient;
typedef FragOrient ChunkOrient;
typedef FragOrient NodeOrient;


/* Fragment positions within a ChunkInstance */
typedef struct {
  CDS_CID_t iid;                // IID of this fragment, used to reference this frag via iidToFragIndex
  CDS_CID_t mateOf;             // the index of the CIFragT of the mate.  Valid even if numLinks > 1
  CDS_CID_t dist;		// index of the DistT record

  CDS_CID_t cid;                // id of the unitig containing this fragment
  CDS_CID_t CIid;               // id of the chunk instance containing this fragment
  LengthT   offset5p;           // offset in containing chunk of suffix
  LengthT   offset3p;           // offset in containing chunk of prefix

  CDS_CID_t contigID;           // id of the containing contig
  LengthT   contigOffset5p;     // offset in containing contig of suffix
  LengthT   contigOffset3p;     // offset in containing contig of prefix

  char      type;               //  a FragType, 5 character values
  char      label;              //  a LabelType, 5 character values

  //  This used to be a union of a struct and an int32; the int32 was
  //  never used, and BPW tried to remove the union....until he
  //  realized that lots of code would need to be changed from
  //  thing.flags.bits.hasMate to thing.hasMate.
  struct {
  struct {
  uint32       hasInternalOnlyCILinks:1;     // If all of this fragment's links are internal to its CI
  uint32       hasInternalOnlyContigLinks:1; // If all of this fragment's links are internal to its Contig
  uint32       isPlaced:1;                   // Fragments is in a contig that is in a scaffold
  uint32       isSingleton:1;                // Singleton unitig
  uint32       isChaff:1;                    // Must be isSingleton, and is not used either directly or indirectly (surrogate) in scaffold
  uint32       innieMate:1;                  // True for regular mate pairs, false for Outtie pairs
  uint32       hasMate:1;                    // True if we have a mate

  LinkType     linkType:8;                   //  mostly unused

  uint32       edgeStatus:7;                 // See edgeStatus field in EdgeCGW_T

  MateStatType mateDetail:8;
  } bits;
  } flags;

  // Enable this to restore compatibility with checkpoints created
  // before Aug 18, 2007.
  //
#if 0
  CDS_CID_t   DEAD_locale;
  SeqInterval DEAD_localePos;
#endif

}CIFragT;

VA_DEF(CIFragT);


static FragOrient getCIFragOrient(CIFragT *frag){
  if(frag->offset3p.mean == frag->offset5p.mean){
    fprintf(stderr,"* Frag %d in unitig %d has 3p=%g 5p=%g\n",frag->iid,frag->cid, frag->offset3p.mean, frag->offset5p.mean);
    assert(0);
  }
  if(frag->offset3p.mean > frag->offset5p.mean)
    return A_B;
  return B_A;
}

static FragOrient GetContigFragOrient(CIFragT *frag){
  if(frag->contigOffset3p.mean == frag->contigOffset5p.mean){
    fprintf(stderr,"* Frag %d in contig %d has 3p=%g 5p=%g\n",frag->iid,frag->contigID, frag->contigOffset3p.mean, frag->contigOffset5p.mean);
    assert(0);
  }
  if(frag->contigOffset3p.mean > frag->contigOffset5p.mean)
    return A_B;
  return B_A;
}


VA_DEF(CDS_COORD_t);

#define CGW_NUM_BUCKETS 3
#define CGW_CUTOFF 5

/* Distance Records */
typedef struct {
  double        mu;              // Calculated from chunk-internal mates
  double        sigma;           // Calculated from chunk-internal mates
  int32         numSamples;      // Redundant -- REMOVE THIS!
  CDS_COORD_t   min;             // Calculated, from contigs
  CDS_COORD_t   max;             // Calculated, from contigs
  int32         bnum;            // number of buckets
  float         bsize;           // size of buckets
  int32        *histogram;
  CDS_COORD_t   lower;
  CDS_COORD_t   upper;
  int32         numReferences;   // Total number of links referencing this distance record
  int32         numBad;
}DistT;

VA_DEF(DistT);

#endif
