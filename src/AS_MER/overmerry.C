
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2006-2007, J. Craig Venter Institute
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

#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
}

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "positionDB.H"

//  The real overlap size:
//  26 - Aiid
//  26 - Biid
//  11 - Apos
//  11 - Bpos
//  10 - k count
//   8 - k size
//   1 - orientation
//   1 - palindrome
//  94 bits, 34 remain


//  The internal overlap
struct kmerhit {
  u64bit   tseq:31;  //  sequence in the table
  u64bit   tpos:11;  //  position in that sequence
  u64bit   qpos:11;  //  position in the query sequence
  u64bit   cnt:10;   //  count of the kmer
  u64bit   ori:1;   //  orientation ; 0 = forward, 1 = reverse
};


int
kmerhitcompare(const void *a, const void *b) {
  const kmerhit *A = (const kmerhit *)a;
  const kmerhit *B = (const kmerhit *)b;
  if (A->tseq < B->tseq)
    return(-1);
  if (A->tseq > B->tseq)
    return(1);
  return(0);
}


inline
u64bit
addHit(chainedSequence *CS, FastASequenceInCore *S, merStream *M,
       kmerhit *&hits, u32bit &hitsLen, u32bit &hitsMax,
       u64bit pos, u64bit cnt,
       u64bit ori) {
  u32bit  seq = CS->sequenceNumberOfPosition(pos);

  pos -= CS->startOf(seq);

  if (S->getIID() != seq) {
    if (hitsLen >= hitsMax) {
      fprintf(stderr, "REALLOC!\n");
    }

    hits[hitsLen].tseq = seq;
    hits[hitsLen].tpos = pos;
    hits[hitsLen].qpos = M->thePositionInSequence();
    hits[hitsLen].cnt  = cnt;
    hits[hitsLen].ori  = ori;
    hitsLen++;

    return(1);
  }
  return(0);
}


int
main(int argc, char **argv) {
  char   *gkpName  = 0L;
  u32bit  merSize  = 23;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      merSize = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkpName == 0L) || (err)) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    exit(1);
  }


  FastABase *gkp = new gkpStoreSequence(gkpName, AS_READ_CLEAR_OBTINI);

  chainedSequence *CS = new chainedSequence;
  CS->setSource(gkp);
  CS->finish();

  merStream    *MS = new merStream(merSize, CS);
  positionDB   *PS = new positionDB(MS, merSize, 0, 26, 0L, 0L, 100, true);

  //  XXXXX: Are we DONE with the MS and CS and gkp here?  Does
  //  positionDB need those still?

  fprintf(stderr, "Go!\n");

  u64bit  *posn    = 0L;
  u64bit   posnMax = 0;
  u64bit   posnLen = 0;

  FastABase            *F = new gkpStoreSequence(gkpName, AS_READ_CLEAR_OBTINI);
  FastASequenceInCore  *S = F->getSequence();

  u64bit  merfound = 0;
  u64bit  ovlfound = 0;

  u32bit     hitsLen = 0;
  u32bit     hitsMax = 1048576;
  kmerhit   *hits    = new kmerhit [hitsMax];

  while (S) {
    merStream  *M = new merStream(merSize, S->sequence(), S->sequenceLength());

    hitsLen = 0;

    while (M->nextMer()) {
      if ((PS->get(M->theFMer(), posn, posnMax, posnLen)) && (posnLen > 1))
        for (u32bit i=0; i<posnLen; i++)
          merfound += addHit(CS, S, M, hits, hitsLen, hitsMax, posn[i], posnLen, 0);
      if ((PS->get(M->theRMer(), posn, posnMax, posnLen)) && (posnLen > 1))
        for (u32bit i=0; i<posnLen; i++)
          merfound += addHit(CS, S, M, hits, hitsLen, hitsMax, posn[i], posnLen, 1);
    }

    //  We have all the hits for this frag.  Sort them by sequence
    //  (the other sequence), then pick out the one with the least
    //  count for each sequence.

    qsort(hits, hitsLen, sizeof(kmerhit), kmerhitcompare);

    u32bit  lowest = 0;

    for (u32bit i=0; i<hitsLen; i++) {
      if ((hits[lowest].tseq != hits[i].tseq) ||
          (i+1 == hitsLen)) {
        ovlfound++;

#if 0
        fprintf(stdout, "seq="u32bitFMT" pos="u64bitFMT" to seq="u64bitFMT" pos="u64bitFMT" count="u64bitFMT" ori=%c\n",
                S->getIID(),
                hits[lowest].qpos,
                hits[lowest].tseq, hits[lowest].tpos,
                hits[lowest].cnt,
                hits[lowest].ori ? 'r' : 'f');
#endif

        lowest = i;
      } else if (hits[lowest].cnt > hits[i].cnt) {
        lowest = i;
      }
    }

    delete M;
    delete S;

    S = F->getSequence();
  }

  fprintf(stderr, "Found "u64bitFMT" mer hits.\n", merfound);
  fprintf(stderr, "Found "u64bitFMT" overlaps.\n", ovlfound);

  delete PS;
  delete MS;
  delete CS;

  delete gkp;
}
