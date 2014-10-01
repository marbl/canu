
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

#ifndef AS_OVS_OVERLAP_H
#define AS_OVS_OVERLAP_H

static const char *rcsid_AS_OVS_OVERLAP_H = "$Id$";

#include "AS_global.H"
#include "AS_MSG_pmesg.H"  //  pretty heavy just to get OverlapMesg.

#define AS_OVS_HNGBITS   (AS_READ_MAX_NORMAL_LEN_BITS + 1)
#define AS_OVS_POSBITS   (AS_READ_MAX_NORMAL_LEN_BITS)
#define AS_OVS_ERRBITS   (12)

#define AS_OVS_OVL_SIZE  (8 + 1 + 2 * AS_OVS_HNGBITS + 2 * AS_OVS_ERRBITS + 2)
#define AS_OVS_MER_SIZE  (3 + 1 + 1 + 2 * AS_OVS_POSBITS + 8 + 8 + 2)
#define AS_OVS_OBT_SIZE  (1 + 4 * AS_OVS_POSBITS + AS_OVS_ERRBITS + 2)

//  Convert q between a condensed/encoded integer and a floating point
//  value.
//
//  Q should be a floating point value between 0.000 and 1.000, and is
//  the fraction error in this alignment.  We are able to encode error
//  up to 0.4000 (40%), with up to four significant figures.
//
//  Previous versions of the overlap store stored any error, with
//  three significant figures, but used 16 bits to do it.  You can get
//  the same effect by using 1000.0 instead of 10000.0 below.

#define AS_OVS_MAX_ERATE          ((1 << AS_OVS_ERRBITS) - 1)

#define AS_OVS_decodeQuality(E)   ((E) / 10000.0)
#define AS_OVS_encodeQuality(Q)   (((Q) < AS_OVS_decodeQuality(AS_OVS_MAX_ERATE)) ? (int)(10000.0 * (Q) + 0.5) : AS_OVS_MAX_ERATE)

#define AS_OVS_TYPE_OVL   0x00
#define AS_OVS_TYPE_OBT   0x01
#define AS_OVS_TYPE_MER   0x02
#define AS_OVS_TYPE_UNS   0x03
#define AS_OVS_TYPE_ANY   0xff

//  These structs are very, very, very carefully laid out.  DO NOT rearrange.
//
//  With 11 or 12 bits in a position, everything fits into 64 bits.  We'd love to use uint32 for
//  this data (to be compatible with the other version of this data), but the OBT data will not fit
//  into 2 32-bit words.
//
//  With more bits, we use a combination of a 32-bit and a 64-bit word.  It's ugly, but having spent
//  an hour trying to layout the data into only 32-bit words, I give up.
//
//  What this means is that everything that uses these data must expect EITHER 32- or 64-bit words.
//  In particular, printf() gets tripped up unless things are cast.
//
//  The "type" field MUST be in the same place in all structs.  We use pad bits to enforce this.
//  Without the padding, the compiler could decide to move things around -- it only need to keep the
//  relative ordering the same.

#if (AS_READ_MAX_NORMAL_LEN_BITS < 13)

#define AS_OVS_NWORDS     2

struct OVSoverlapOVL {
  uint64  datpad1            :64 - 1 - 2 * AS_OVS_HNGBITS - 2 * AS_OVS_ERRBITS - 8 - 2;
  uint64  flipped            :1;
  int64   a_hang             :AS_OVS_HNGBITS;
  int64   b_hang             :AS_OVS_HNGBITS;
  uint64  orig_erate         :AS_OVS_ERRBITS;
  uint64  corr_erate         :AS_OVS_ERRBITS;
  uint64  seed_value         :8;
  uint64  type               :2;
};

struct OVSoverlapMER {
  uint64  datpad1            :64 - 1 - 1 - 2 * AS_OVS_POSBITS - 3 - 8 - 8 - 2;
  uint64  fwd                :1;
  uint64  palindrome         :1;
  uint64  a_pos              :AS_OVS_POSBITS;
  uint64  b_pos              :AS_OVS_POSBITS;
  uint64  compression_length :3;
  uint64  k_count            :8;
  uint64  k_len              :8;
  uint64  type               :2;
};

struct OVSoverlapOBT {
  uint64  datpad1            :64 - 1 - 4 * AS_OVS_POSBITS - AS_OVS_ERRBITS - 2;
  uint64  fwd                :1;
  uint64  a_beg              :AS_OVS_POSBITS;
  uint64  a_end              :AS_OVS_POSBITS;
  uint64  b_beg              :AS_OVS_POSBITS;
  uint64  b_end_hi           :AS_OVS_POSBITS - 9;
  uint64  b_end_lo           :9;
  uint64  erate              :AS_OVS_ERRBITS;
  uint64  type               :2;
};

#elif ((AS_READ_MAX_NORMAL_LEN_BITS == 13) && defined(OVS_DANGEROUSLY_OVERSIZE))
#warning OVSE_DANGEROUSLY_OVERSIZE defined; OBT will not work

//  In this mode, we can max out the space for OVSoverlapOVL, but OVSoverlapOBT is not supported.
//  To make it fit in 2 words, we could remove one bit from each of the position fields, or 4 bits
//  from the erate.  But to make it completely obvious this is broken, we instead reduce the
//  position fields to exactly one bit.  OBT will still (probably) run, but reads will be
//  obliterated.

#define AS_OVS_NWORDS     2

struct OVSoverlapOVL {
  uint64  datpad1            :64 - 1 - 2 * AS_OVS_HNGBITS - 2 * AS_OVS_ERRBITS - 8 - 2;
  uint64  flipped            :1;
  int64   a_hang             :AS_OVS_HNGBITS;
  int64   b_hang             :AS_OVS_HNGBITS;
  uint64  orig_erate         :AS_OVS_ERRBITS;
  uint64  corr_erate         :AS_OVS_ERRBITS;
  uint64  seed_value         :8;
  uint64  type               :2;
};

struct OVSoverlapMER {
  uint64  datpad1            :64 - 1 - 1 - 2 * AS_OVS_POSBITS - 3 - 8 - 8 - 2;
  uint64  fwd                :1;
  uint64  palindrome         :1;
  uint64  a_pos              :AS_OVS_POSBITS;
  uint64  b_pos              :AS_OVS_POSBITS;
  uint64  compression_length :3;
  uint64  k_count            :8;
  uint64  k_len              :8;
  uint64  type               :2;
};

struct OVSoverlapOBT {
  uint64  datpad1            :64 - 1 - 4 * (1) - 9 - AS_OVS_ERRBITS - 2;
  uint64  fwd                :1;
  uint64  a_beg              :1;
  uint64  a_end              :1;
  uint64  b_beg              :1;
  uint64  b_end_hi           :1;
  uint64  b_end_lo           :9;
  uint64  erate              :AS_OVS_ERRBITS;
  uint64  type               :2;
};

#else

#define AS_OVS_NWORDS 3

struct OVSoverlapOVL {
  uint32  seed_value         :8;
  uint32  orig_erate         :AS_OVS_ERRBITS;
  uint32  corr_erate         :AS_OVS_ERRBITS;

  uint32  datpad1            :32 - AS_OVS_HNGBITS;
  int32   a_hang             :AS_OVS_HNGBITS;

  uint32  datpad2            :32 - AS_OVS_HNGBITS - 1 - 2;
  int32   b_hang             :AS_OVS_HNGBITS;
  uint32  flipped            :1;
  uint32  type               :2;
};

struct OVSoverlapMER {
  uint32  datpad1            :32 - AS_OVS_POSBITS;
  uint32  a_pos              :AS_OVS_POSBITS;

  uint32  datpad2            :32 - AS_OVS_POSBITS;
  uint32  b_pos              :AS_OVS_POSBITS;

  uint32  datpad3            :32 - 3 - 8 - 8 - 1 - 1 - 2;
  uint32  compression_length :3;
  uint32  k_count            :8;
  uint32  k_len              :8;
  uint32  palindrome         :1;
  uint32  fwd                :1;
  uint32  type               :2;
};

struct OVSoverlapOBT {
#if (32 - AS_OVS_ERRBITS - AS_OVS_POSBITS > 0)
  uint32  datpad1            :32 - AS_OVS_ERRBITS - AS_OVS_POSBITS;
#endif
  uint32  erate              :AS_OVS_ERRBITS;
  uint32  a_beg              :AS_OVS_POSBITS;

  uint32  datpad2            :32 - 2 * AS_OVS_POSBITS + 9;
  uint32  a_end              :AS_OVS_POSBITS;
  uint32  b_end_hi           :AS_OVS_POSBITS - 9;

#if (32 - AS_OVS_POSBITS - 9 - 1 - 2 > 0)
  uint32  datpad3            :32 - AS_OVS_POSBITS - 9 - 1 - 2;
#endif
  uint32  b_beg              :AS_OVS_POSBITS;
  uint32  b_end_lo           :9;
  uint32  fwd                :1;
  uint32  type               :2;
};

#endif

#define AS_OVS_WORD_SIZE sizeof(uint32)

typedef union {
  uint32                dat[AS_OVS_NWORDS];
  struct OVSoverlapOVL  ovl;
  struct OVSoverlapMER  mer;
  struct OVSoverlapOBT  obt;
} OVSoverlapDAT;

typedef struct {
  AS_IID         b_iid;
  OVSoverlapDAT  dat;
} OVSoverlapINT;

typedef struct {
  AS_IID         a_iid;
  AS_IID         b_iid;
  OVSoverlapDAT  dat;
} OVSoverlap;

inline
bool
operator<(const OVSoverlap &A, const OVSoverlap &B) {
  if (A.a_iid      < B.a_iid)       return(true);
  if (A.a_iid      > B.a_iid)       return(false);
  if (A.b_iid      < B.b_iid)       return(true);
  if (A.b_iid      > B.b_iid)       return(false);
  if (A.dat.dat[0] < B.dat.dat[0])  return(true);
  if (A.dat.dat[0] > B.dat.dat[0])  return(false);
  if (A.dat.dat[1] < B.dat.dat[1])  return(true);
  if (A.dat.dat[1] > B.dat.dat[1])  return(false);
#if (AS_OVS_NWORDS > 2)
  if (A.dat.dat[2] < B.dat.dat[2])  return(true);
  if (A.dat.dat[2] > B.dat.dat[2])  return(false);
#endif
#if (AS_OVS_NWORDS > 3)
  if (A.dat.dat[3] < B.dat.dat[3])  return(true);
  if (A.dat.dat[3] > B.dat.dat[3])  return(false);
#endif
  return(false);
}


void  AS_OVS_convertOverlapMesgToOVSoverlap(OverlapMesg *omesg, OVSoverlap *ovs);
int   AS_OVS_convertOVLdumpToOVSoverlap(char *line, OVSoverlap *olap);
int   AS_OVS_convertOBTdumpToOVSoverlap(char *line, OVSoverlap *olap);

void  AS_OVS_convertOVLoverlapToOBToverlap(OVSoverlap &overlap, uint32 clrLenA, uint32 clrLenB);

char *AS_OVS_toString(char *outstr, OVSoverlap &olap);


static
uint32
AS_OVS_overlapAEndIs5prime(const OVSoverlap& olap) {
  return((olap.dat.ovl.a_hang < 0) && (olap.dat.ovl.b_hang < 0));
}

static
uint32
AS_OVS_overlapAEndIs3prime(const OVSoverlap& olap) {
  return((olap.dat.ovl.a_hang > 0) && (olap.dat.ovl.b_hang > 0));
}

static
uint32
AS_OVS_overlapAIsContained(const OVSoverlap& olap) {
  return((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0));
}

static
uint32
AS_OVS_overlapAIsContainer(const OVSoverlap& olap) {
  return((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0));
}



static
uint32
AS_OVS_overlapBEndIs5prime(const OVSoverlap& olap) {
  return((AS_OVS_overlapAEndIs5prime(olap) && (olap.dat.ovl.flipped == true)) ||
         (AS_OVS_overlapAEndIs3prime(olap) && (olap.dat.ovl.flipped == false)));
}

static
uint32
AS_OVS_overlapBEndIs3prime(const OVSoverlap& olap) {
  return((AS_OVS_overlapAEndIs5prime(olap) && (olap.dat.ovl.flipped == false)) ||
         (AS_OVS_overlapAEndIs3prime(olap) && (olap.dat.ovl.flipped == true)));
}

static
uint32
AS_OVS_overlapBIsContained(const OVSoverlap& olap) {
  return((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0));
}

static
uint32
AS_OVS_overlapBIsContainer(const OVSoverlap& olap) {
  return((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0));
}


static
void
AS_OVS_testEndComparisons_Test(const OVSoverlap &olap, uint32 c=1) {

  if (AS_OVS_overlapAEndIs5prime(olap) + AS_OVS_overlapAEndIs3prime(olap) + AS_OVS_overlapAIsContained(olap) + AS_OVS_overlapAIsContainer(olap) != c) {
    fprintf(stderr, "A AS_OVS_overlapAEndIs5prime(olap) = %d\n", AS_OVS_overlapAEndIs5prime(olap));
    fprintf(stderr, "A AS_OVS_overlapAEndIs3prime(olap) = %d\n", AS_OVS_overlapAEndIs3prime(olap));
    fprintf(stderr, "A AS_OVS_overlapAIsContained(olap) = %d\n", AS_OVS_overlapAIsContained(olap));
    fprintf(stderr, "A AS_OVS_overlapAIsContainer(olap) = %d\n", AS_OVS_overlapAIsContainer(olap));
    fprintf(stderr, "a_hang "F_S64" b_hang "F_S64" flipped "F_U64"\n", (int64)olap.dat.ovl.a_hang, (int64)olap.dat.ovl.b_hang, (int64)olap.dat.ovl.flipped);
  }

  if (AS_OVS_overlapBEndIs5prime(olap) + AS_OVS_overlapBEndIs3prime(olap) + AS_OVS_overlapBIsContained(olap) + AS_OVS_overlapBIsContainer(olap) != c) {
    fprintf(stderr, "B AS_OVS_overlapBEndIs5prime(olap) = %d\n", AS_OVS_overlapBEndIs5prime(olap));
    fprintf(stderr, "B AS_OVS_overlapBEndIs3prime(olap) = %d\n", AS_OVS_overlapBEndIs3prime(olap));
    fprintf(stderr, "B AS_OVS_overlapBIsContained(olap) = %d\n", AS_OVS_overlapBIsContained(olap));
    fprintf(stderr, "B AS_OVS_overlapBIsContainer(olap) = %d\n", AS_OVS_overlapBIsContainer(olap));
    fprintf(stderr, "a_hang "F_S64" b_hang "F_S64" flipped "F_U64"\n", (int64)olap.dat.ovl.a_hang, (int64)olap.dat.ovl.b_hang, (int64)olap.dat.ovl.flipped);
  }

  assert(AS_OVS_overlapAEndIs5prime(olap) + AS_OVS_overlapAEndIs3prime(olap) + AS_OVS_overlapAIsContained(olap) + AS_OVS_overlapAIsContainer(olap) == c);
  assert(AS_OVS_overlapBEndIs5prime(olap) + AS_OVS_overlapBEndIs3prime(olap) + AS_OVS_overlapBIsContained(olap) + AS_OVS_overlapBIsContainer(olap) == c);
}

static
void
AS_OVS_testEndComparisons(void) {
  OVSoverlap  olap;

  olap.dat.ovl.a_hang = -1;  olap.dat.ovl.b_hang = -1;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang = -1;  olap.dat.ovl.b_hang =  0;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang = -1;  olap.dat.ovl.b_hang =  1;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);

  olap.dat.ovl.a_hang =  0;  olap.dat.ovl.b_hang = -1;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang =  0;  olap.dat.ovl.b_hang =  0;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap, 2);
  olap.dat.ovl.a_hang =  0;  olap.dat.ovl.b_hang =  1;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);

  olap.dat.ovl.a_hang =  1;  olap.dat.ovl.b_hang = -1;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang =  1;  olap.dat.ovl.b_hang =  0;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang =  1;  olap.dat.ovl.b_hang =  1;  olap.dat.ovl.flipped = 0;  AS_OVS_testEndComparisons_Test(olap);


  olap.dat.ovl.a_hang = -1;  olap.dat.ovl.b_hang = -1;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang = -1;  olap.dat.ovl.b_hang =  0;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang = -1;  olap.dat.ovl.b_hang =  1;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);

  olap.dat.ovl.a_hang =  0;  olap.dat.ovl.b_hang = -1;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang =  0;  olap.dat.ovl.b_hang =  0;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap, 2);
  olap.dat.ovl.a_hang =  0;  olap.dat.ovl.b_hang =  1;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);

  olap.dat.ovl.a_hang =  1;  olap.dat.ovl.b_hang = -1;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang =  1;  olap.dat.ovl.b_hang =  0;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);
  olap.dat.ovl.a_hang =  1;  olap.dat.ovl.b_hang =  1;  olap.dat.ovl.flipped = 1;  AS_OVS_testEndComparisons_Test(olap);
}



#endif  //  AS_OVS_OVERLAP_H
