
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
#include "AS_OVS_overlapStore.h"
}

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "positionDB.H"


//  Instead of using The internal overlap, which has enough extra
//  stuff in it that we cannot store a sequence iid for the table
//  sequence, we need to make our own overlap structure.
//
struct kmerhit {
  u64bit   tseq:30;              //  sequence in the table
  u64bit   tpos:AS_OVS_POSBITS;  //  position in that sequence
  u64bit   qpos:AS_OVS_POSBITS;  //  position in the query sequence
  u64bit   cnt:8;                //  count of the kmer
  u64bit   pal:1;                //  palindromic ; 0 = nope,    1 = yup
  u64bit   fwd:1;                //  orientation ; 0 = reverse, 1 = forward
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
addHit(chainedSequence *CS, seqInCore *S, merStream *M,
       kmerhit *&hits, u32bit &hitsLen, u32bit &hitsMax,
       u64bit pos, u64bit cnt,
       u64bit pal,
       u64bit fwd) {
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
    hits[hitsLen].pal  = pal;
    hits[hitsLen].fwd  = fwd;
    hitsLen++;

    return(1);
  }
  return(0);
}


int
main(int argc, char **argv) {
  char   *gkpPath  = 0L;
  char    gkpName[FILENAME_MAX + 64] = {0};
  u32bit  merSize  = 23;

  assert(sizeof(kmerhit) == 8);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpPath = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      merSize = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkpPath == 0L) || (err)) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    exit(1);
  }


  seqFactory::instance()->registerFile(new gkpStoreSequence());

  //  Open the gatekeeper store as a kmer seqFile.  We need to
  //  dynamic_cast this back to our gkpStoreSequence, so we can access
  //  methods defined only on that object.
  //
  sprintf(gkpName, "%s:obtini", gkpPath);
  gkpStoreSequence *gkpseq = dynamic_cast<gkpStoreSequence*>(openSeqFile(gkpName));
  if (gkpseq == 0L) {
    fprintf(stderr, "%s: invalid input file '%s' (not a GateKeeperStore?).\n", gkpName);
    exit(1);
  }


  //  XXXX  extend merstream with spaced, compressed, skips, transitions


#if 1
  u32bit i;
  seqOnDisk *s;

  i=0;
  gkpseq->find(i);
  s = gkpseq->getSequenceOnDisk();
  fprintf(stderr, "0: '%s'\n", s->sequence());

  i=1;
  gkpseq->find(i);
  s = gkpseq->getSequenceOnDisk();
  fprintf(stderr, "1: '%s'\n", s->sequence());
#endif

  chainedSequence *CS = new chainedSequence;
  CS->setSource(gkpseq);
  CS->finish();

  merStream    *MS = new merStream(merSize, CS);
  positionDB   *PS = new positionDB(MS, merSize, 0, 26, 0L, 0L, 100, true);

  u64bit  *posn    = 0L;
  u64bit   posnMax = 0;
  u64bit   posnLen = 0;

  fragRecord       *fr  = new_fragRecord();
  GateKeeperStore  *gkp = openGateKeeperStore(gkpPath, FALSE);
  FragStream       *frg = openFragStream(gkp, FRAG_S_SEQ);

  u64bit  merfound = 0;
  u64bit  ovlfound = 0;

  u32bit     hitsLen = 0;
  u32bit     hitsMax = 1048576;
  kmerhit   *hits    = new kmerhit [hitsMax];

  fprintf(stderr, "Go!\n");

#ifdef BINARYOUTPUY
  BinaryOverlapFile *binout = AS_OVS_createBinaryOverlapFile("-", FALSE);
  OVSoverlap         overlap;
#endif

  while (nextFragStream(frg, fr)) {
    char        *seq = new char [getFragRecordSequenceLength(fr) + 1];
    uint32       beg = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBTINI);
    uint32       end = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBTINI);
    uint32       len = end - beg;

    strncpy(seq, getFragRecordSequence(fr) + beg, len);
    seq[len] = 0;

    seqInCore   *S   = new seqInCore(getFragRecordIID(fr), 0L, 0, seq, len);
    merStream  *M    = new merStream(merSize, S);

    hitsLen = 0;

    while (M->nextMer()) {
      if (M->theFMer() == M->theRMer()) {
        if ((PS->get(M->theFMer(), posn, posnMax, posnLen)) && (posnLen > 1))
          for (u32bit i=0; i<posnLen; i++)
            merfound += addHit(CS, S, M, hits, hitsLen, hitsMax, posn[i], posnLen, 1, 0);
      } else {
        if ((PS->get(M->theFMer(), posn, posnMax, posnLen)) && (posnLen > 1))
          for (u32bit i=0; i<posnLen; i++)
            merfound += addHit(CS, S, M, hits, hitsLen, hitsMax, posn[i], posnLen, 0, 1);
        if ((PS->get(M->theRMer(), posn, posnMax, posnLen)) && (posnLen > 1))
          for (u32bit i=0; i<posnLen; i++)
            merfound += addHit(CS, S, M, hits, hitsLen, hitsMax, posn[i], posnLen, 0, 0);
      }
    }

    //  We have all the hits for this frag.  Sort them by sequence
    //  (the other sequence), then pick out the one with the least
    //  count for each sequence.

    qsort(hits, hitsLen, sizeof(kmerhit), kmerhitcompare);

    u32bit  lowest = 0;

    for (u32bit i=0; i<=hitsLen; i++) {

#if 0
      //  Debug, I guess.  Generates lots of output, since frags with
      //  big identical overlaps will have lots and lots of mers in
      //  common.
      //
      if (i != hitsLen) {
        fprintf(stdout, u32bitFMT"\t"u64bitFMT"\t"u32bitFMT"\t"u64bitFMT"\t%c\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\tTAG\n",
                S->getIID(), hits[i].qpos,
                hits[i].tseq, hits[i].tpos,
                hits[i].pal ? 'p' : (hits[i].fwd ? 'f' : 'r'),
                0,
                hits[i].cnt,
                merSize);
      }
#endif

      if ((i == hitsLen) || (hits[lowest].tseq != hits[i].tseq)) {

        //  We either found a different sequence than the one we were
        //  looking at, or we've gone past the end.

        ovlfound++;

        if (hits[lowest].fwd == 0)
          hits[lowest].qpos = S->sequenceLength() - hits[lowest].qpos;

        //  Output IID1, pos1, IID2, pos2, orient/palindrome, compressionLength, count, kmersize

        //  Adjust coords to be relative to whole read
        //
        hits[lowest].tpos += gkpseq->clrBeg(hits[lowest].tseq);
        hits[lowest].qpos += gkpseq->clrBeg(S->getIID());

        //  And reverse if needed.
        //
        if (hits[lowest].fwd == false)
          hits[lowest].qpos = gkpseq->sequenceLength(S->getIID()) - hits[lowest].qpos - merSize;

#ifdef BINARYOUTPUT
        overlap.a_iid              = hits[lowest].tseq;
        overlap.b_iid              = S->getIID();
        overlap.datpad             = 0;
        overlap.compression_length = 0;
        overlap.fwd                = hits[lowest].fwd;
        overlap.palindrome         = hits[lowest].pal;
        overlap.a_pos              = hits[lowest].tpos;
        overlap.b_pos              = hits[lowest].qpos;
        overlap.k_count            = hits[lowest].cnt;
        overlap.k_len              = merSize;
        overlap.type               = AS_OVS_TYPE_MER;
        
        AS_OVS_writeOverlap(binout, &overlap);
#else
        fprintf(stdout, "%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\n",
                (int)hits[lowest].tseq,
                (int)hits[lowest].tpos,
                (int)S->getIID(),
                (int)hits[lowest].qpos,
                hits[lowest].pal ? 'p' : (hits[lowest].fwd ? 'f' : 'r'),
                (int)0,
                (int)hits[lowest].cnt,
                (int)merSize);
#endif

        lowest = i;
      } else if (hits[lowest].cnt > hits[i].cnt) {
        lowest = i;
      }
    }

    delete M;
    delete S;
  }

  fprintf(stderr, "Found "u64bitFMT" mer hits.\n", merfound);
  fprintf(stderr, "Found "u64bitFMT" overlaps.\n", ovlfound);

#ifdef BINARYOUTPUT
  AS_OVS_closeBinaryOverlapFile(binout);
#endif

  closeFragStream(frg);
  closeGateKeeperStore(gkp);

  delete PS;
  delete MS;
  delete CS;

  delete gkpseq;
}
