
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

#ifndef INPUTDATATYPES_CGW_H
#define INPUTDATATYPES_CGW_H

static const char *rcsid_INPUTDATATYPES_CGW_H = "$Id: InputDataTypes_CGW.h,v 1.26 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

typedef struct {
  double  mean;
  double variance;
}LengthT;


typedef enum {
  AS_SINGLETON    = (int)'S',
  AS_INTERCHUNK_A = (int)'A',
  AS_INTERCHUNK_B = (int)'B',
  AS_INTRACHUNK   = (int)'I'
} LabelType;


/* Fragment positions within a ChunkInstance */
typedef struct {
  CDS_CID_t read_iid;           // IID of this fragment
  CDS_CID_t mate_iid;           // IID of the mate
  CDS_CID_t dist;		// index of the DistT record

  CDS_CID_t cid;                // id of the unitig containing this fragment
  CDS_CID_t CIid;               // id of the chunk instance containing this fragment
  CDS_CID_t contigID;           // id of the containing contig

  LengthT   offset5p;           // offset in containing chunk of suffix
  LengthT   offset3p;           // offset in containing chunk of prefix

  LengthT   contigOffset5p;     // offset in containing contig of suffix
  LengthT   contigOffset3p;     // offset in containing contig of prefix

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
  uint32       isDeleted:1;                  // Frag exists in store, but marked deleted
  uint32       innieMate:1;                  // True for regular mate pairs, false for Outtie pairs
  uint32       hasMate:1;                    // True if we have a mate

  uint32       edgeStatus:7;                 // See edgeStatus field in EdgeCGW_T
  uint32       chunkLabel:8;                 // A LabelType

  MateStatType mateDetail:8;
  } bits;
  } flags;
}CIFragT;

VA_DEF(CIFragT)


static SequenceOrient getCIFragOrient(CIFragT *frag){
  SequenceOrient  orient;

  if(frag->offset3p.mean == frag->offset5p.mean){
    fprintf(stderr,"* Frag %d in unitig %d has 3p=%g 5p=%g\n",frag->read_iid,frag->cid, frag->offset3p.mean, frag->offset5p.mean);
    assert(0);
  }
  orient.setIsForward(frag->offset5p.mean < frag->offset3p.mean);
  return(orient);
}

static SequenceOrient GetContigFragOrient(CIFragT *frag){
  SequenceOrient  orient;

  if(frag->contigOffset3p.mean == frag->contigOffset5p.mean){
    fprintf(stderr,"* Frag %d in contig %d has 3p=%g 5p=%g\n",frag->read_iid,frag->contigID, frag->contigOffset3p.mean, frag->contigOffset5p.mean);
    assert(0);
  }
  orient.setIsForward(frag->contigOffset5p.mean < frag->contigOffset3p.mean);
  return(orient);
}

#define CGW_NUM_BUCKETS 3
#define CGW_CUTOFF 5

/* Distance Records */
typedef struct {
  double        mu;              // Calculated from chunk-internal mates
  double        sigma;           // Calculated from chunk-internal mates
  int32         numSamples;      // Redundant -- REMOVE THIS!
  int32         min;             // Calculated, from contigs
  int32         max;             // Calculated, from contigs
  int32         bnum;            // number of buckets
  double        bsize;           // size of buckets
  int32        *histogram;
  int32         lower;
  int32         upper;
  int32         numReferences;   // Total number of links referencing this distance record
  int32         numBad;
}DistT;

VA_DEF(DistT)

#endif
