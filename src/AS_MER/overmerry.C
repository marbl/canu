
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


//  using external mer counts - we'll need to extend the posDB to have
//  more stuff in it, or construct another mer lookup hash.
//
//  parallel searching
//
//  subset of gkpstore
//
//  extend merstream with spaced, compressed, skips, transitions
//
//  M039 takes 134u total, 30u to build the table.


#undef EXPENSIVE_CHECK
#ifdef EXPENSIVE_CHECK

//  For debug, swiped from AS_OVL
/*  Return the DNA complement of  ch . */
static char  Complement (char ch) {
   switch  (((int) ch)) {
      case  'A' : return  'T';
      case  'C' : return  'G';
      case  'G' : return  'C';
      case  'T' : return  'A';
      case  'N' : return  'N';
      default :
        fprintf (stderr, "ERROR(complement):  Unexpected character `%c\'\n", ch);
        exit (-1);
     }
   return  'X';    // Just to make the compiler happy
  }


/* Set string  s  to its DNA reverse complement. */
static void  Rev_Complement (char * s) {
   char  ch;
   int  i, j, len; 
   len = strlen (s); 
   for  (i = 0, j = len - 1;  i < j;  i ++, j --) {
      ch = Complement (s [i]);
      s [i] = Complement (s [j]);
      s [j] = ch;
     } 
   if  (i == j) s [i] = Complement (s [i]); 
  }

#endif





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
  u64bit   pad:2;
};



//  Sort by increasing iid, and increasing count.
int
kmerhitcompare(const void *a, const void *b) {
  const kmerhit *A = (const kmerhit *)a;
  const kmerhit *B = (const kmerhit *)b;
  if (A->tseq < B->tseq)
    return(-1);
  if (A->tseq > B->tseq)
    return(1);
  if (A->cnt < B->cnt)
    return(-1);
  if (A->cnt > B->cnt)
    return(1);
  return(0);
}


inline
u64bit
addHit(chainedSequence *CS, CDS_IID_t iid, merStream *M,
       kmerhit *&hits, u32bit &hitsLen, u32bit &hitsMax,
       u64bit pos, u64bit cnt,
       u64bit pal,
       u64bit fwd) {
  u32bit  seq = CS->sequenceNumberOfPosition(pos);

  pos -= CS->startOf(seq);
  seq  = CS->IIDOf(seq);

  if (iid == seq)
    return(0);

  if (iid < seq)
    return(0);

  if (hitsLen >= hitsMax) {
    if (hitsMax == 0) {
      hitsMax = 256;
      hits    = new kmerhit [hitsMax];
    } else {
      hitsMax *= 2;
      kmerhit *h = new kmerhit [hitsMax];
      memcpy(h, hits, sizeof(kmerhit) * hitsLen);
      delete [] hits;
      hits = h;
      fprintf(stderr, "reallocating to hitsMax "u32bitFMT"\n", hitsMax);
    }
  }

  hits[hitsLen].tseq = seq;
  hits[hitsLen].tpos = pos;
  hits[hitsLen].qpos = M->thePositionInSequence();
  hits[hitsLen].cnt  = cnt;
  hits[hitsLen].pal  = pal;
  hits[hitsLen].fwd  = fwd;

#if 0
  if (((seq == 10) && (iid ==  1)) ||
      ((seq ==  1) && (iid == 10)))
    fprintf(stderr, "addHit: "u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t%c\n",
            hits[hitsLen].tseq,
            hits[hitsLen].tpos,
            hits[hitsLen].qpos,
            hits[hitsLen].cnt,
            hits[hitsLen].pal ? 'p' : (hits[hitsLen].fwd ? 'f' : 'r'));
#endif

  hitsLen++;

  return(1);
}


int
main(int argc, char **argv) {
  char   *gkpPath    = 0L;
  char    gkpName[FILENAME_MAX + 64] = {0};
  u32bit  merSize    = 23;
  u32bit  maxCount   = 40;
  char   *outputName = NULL;

  assert(sizeof(kmerhit) == 8);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpPath = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      merSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-o") == 0) {
      outputName = argv[++arg];
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkpPath == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-g gkpStore] [-m merSize] [-o outputName]\n", argv[0]);
    exit(1);
  }

  seqFactory::instance()->registerFile(new gkpStoreSequence());

  //  Open the gatekeeper store as a kmer seqFile.  We need to
  //  dynamic_cast this back to our gkpStoreSequence, so we can access
  //  methods defined only on that object.
  //
  sprintf(gkpName, "%s:obt", gkpPath);
  gkpStoreSequence *gkpseq = dynamic_cast<gkpStoreSequence*>(openSeqFile(gkpName));
  if (gkpseq == 0L) {
    fprintf(stderr, "%s: invalid input file '%s' (not a GateKeeperStore?).\n", gkpName);
    exit(1);
  }

  chainedSequence *CS = new chainedSequence;
  CS->setSeparatorLength(5);
  CS->setSource(gkpseq);
  CS->finish();

  merStream    *MS = new merStream(merSize, CS);
  positionDB   *PS = new positionDB(MS, merSize, 0, 20, 0L, 0L, maxCount, true);

  u64bit  *posnF    = 0L;
  u64bit   posnFMax = 0;
  u64bit   posnFLen = 0;

  u64bit  *posnR    = 0L;
  u64bit   posnRMax = 0;
  u64bit   posnRLen = 0;

  u64bit  merfound = 0;
  u64bit  ovlfound = 0;

  u32bit     hitsLen = 0;
  u32bit     hitsMax = 0;
  kmerhit   *hits    = 0;

  BinaryOverlapFile *binout = AS_OVS_createBinaryOverlapFile(outputName, FALSE);
  OVSoverlap         overlap;

#ifdef EXPENSIVE_CHECK
  GateKeeperStore  *dbggkp = openGateKeeperStore(gkpPath, FALSE);
  fragRecord       *dbgf1  = new_fragRecord();
  fragRecord       *dbgf2  = new_fragRecord();
  char              dbgs1[AS_READ_MAX_LEN+1];
  char              dbgs2[AS_READ_MAX_LEN+1];
#endif

  fragRecord       *fr  = new_fragRecord();
  GateKeeperStore  *gkp = openGateKeeperStore(gkpPath, FALSE);
  FragStream       *frs = openFragStream(gkp, FRAG_S_SEQ);


  while (nextFragStream(frs, fr)) {
    fprintf(stderr, "WORKING on "F_UID","F_IID"\r",
            getFragRecordUID(fr), getFragRecordIID(fr));

    if (getFragRecordIsDeleted(fr))
      continue;

    uint32       beg = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBT);
    uint32       end = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBT);
    uint32       tln = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_UNTRIM);
    CDS_IID_t    iid = getFragRecordIID(fr);
    merStream   *M   = new merStream(merSize, getFragRecordSequence(fr), beg, end - beg);

    hitsLen = 0;

    //  We can "simulate" a canonical mercount by getting the count
    //  for both forward and reverse.

    while (M->nextMer()) {
      if (M->theFMer() == M->theRMer()) {
        PS->get(M->theFMer(), posnF, posnFMax, posnFLen);
        if (posnFLen > 1)
          for (u32bit i=0; i<posnFLen; i++)
            merfound += addHit(CS, getFragRecordIID(fr), M, hits, hitsLen, hitsMax, posnF[i], posnFLen, 1, 0);
      } else {
        PS->get(M->theFMer(), posnF, posnFMax, posnFLen);
        PS->get(M->theRMer(), posnR, posnRMax, posnRLen);

        u64bit totalLen = posnFLen + posnRLen - 1;  //  the canonical mer count

        if (posnFLen > 1)
          for (u32bit i=0; i<posnFLen; i++)
            merfound += addHit(CS, iid, M, hits, hitsLen, hitsMax, posnF[i], totalLen, 0, 1);
        if (posnRLen > 1)
          for (u32bit i=0; i<posnRLen; i++)
            merfound += addHit(CS, iid, M, hits, hitsLen, hitsMax, posnR[i], totalLen, 0, 0);
      }
    }

    if (hitsLen == 0)
      continue;

    //  We have all the hits for this frag.  Sort them by sequence
    //  (the other sequence), then pick out the one with the least
    //  count for each sequence.

    qsort(hits, hitsLen, sizeof(kmerhit), kmerhitcompare);

#if 0
    //  Debug, I guess.  Generates lots of output, since frags with
    //  big identical overlaps will have lots and lots of mers in
    //  common.
    //
    for (u32bit i=0; i<hitsLen; i++) {
      if (i != hitsLen) {
        fprintf(stderr, u32bitFMT"\t"u64bitFMT"\t"u32bitFMT"\t"u64bitFMT"\t%c\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\tTAG\n",
                hits[i].tseq, hits[i].tpos,
                iid,  hits[i].qpos,
                hits[i].pal ? 'p' : (hits[i].fwd ? 'f' : 'r'),
                0,
                hits[i].cnt,
                merSize);
      }
    }
#endif



    for (u32bit i=0; i<hitsLen; i++) {

      //  By the definition of our sort, the least common mer is the
      //  first hit in the list for each pair of sequences.
      //
      ovlfound++;

      //  Adjust coords to be relative to whole read
      //
      hits[i].tpos += gkpseq->clrBeg(hits[i].tseq);
      hits[i].qpos += beg;

      //  Reverse if needed
      //
      if (hits[i].fwd == false)
        hits[i].qpos = tln - hits[i].qpos - merSize;

      overlap.a_iid                      = hits[i].tseq;
      overlap.b_iid                      = iid;
      overlap.dat.mer.datpad             = 0;
      overlap.dat.mer.compression_length = 0;
      overlap.dat.mer.fwd                = hits[i].fwd;
      overlap.dat.mer.palindrome         = hits[i].pal;
      overlap.dat.mer.a_pos              = hits[i].tpos;
      overlap.dat.mer.b_pos              = hits[i].qpos;
      overlap.dat.mer.k_count            = hits[i].cnt;
      overlap.dat.mer.k_len              = merSize;
      overlap.dat.mer.type               = AS_OVS_TYPE_MER;

#if 0
      //  some ascii output
      fprintf(stderr, "%d\t%d\t%c\t%d\t%d\n",
              (int)overlap.a_iid,
              (int)overlap.dat.mer.a_pos,
              overlap.dat.mer.palindrome ? 'p' : (overlap.dat.mer.fwd ? 'f' : 'r'),
              (int)overlap.b_iid,
              (int)overlap.dat.mer.b_pos);
#endif

      AS_OVS_writeOverlap(binout, &overlap);


#ifdef EXPENSIVE_CHECK
      getFrag(dbggkp, overlap.a_iid, dbgf1, FRAG_S_SEQ);
      getFrag(dbggkp, overlap.b_iid, dbgf2, FRAG_S_SEQ);

      //  Debugging....lots of it....
      //fprintf(stderr, ">"F_IID"\n%s\n", getFragRecordIID(dbgf1), getFragRecordSequence(dbgf1));
      //fprintf(stderr, ">"F_IID"\n%s\n", getFragRecordIID(dbgf2), getFragRecordSequence(dbgf2));

      if (overlap.dat.mer.fwd == false)
        Rev_Complement(getFragRecordSequence(dbgf2));

      if (strncasecmp(getFragRecordSequence(dbgf1) + overlap.dat.mer.a_pos,
                      getFragRecordSequence(dbgf2) + overlap.dat.mer.b_pos, merSize) != 0) {

        fprintf(stderr, ">"F_IID"\n%s\n", getFragRecordIID(dbgf1), getFragRecordSequence(dbgf1));
        fprintf(stderr, ">"F_IID"\n%s\n", getFragRecordIID(dbgf2), getFragRecordSequence(dbgf2));

        fprintf(stderr, "%d\t%d\t%c\t%d\t%d\n",
                (int)overlap.a_iid,
                (int)overlap.dat.mer.a_pos,
                overlap.dat.mer.palindrome ? 'p' : (overlap.dat.mer.fwd ? 'f' : 'r'),
                (int)overlap.b_iid,
                (int)overlap.dat.mer.b_pos);

        char  *a = getFragRecordSequence(dbgf1) + overlap.dat.mer.a_pos;
        char  *b = getFragRecordSequence(dbgf2);
        int    i = 0;
        for (; *b; b++, i++)
          if (strncasecmp(a, b, overlap.dat.mer.k_len) == 0)
            fprintf(stderr, "found at %d\n", i);

        getFragRecordSequence(dbgf1)[overlap.dat.mer.a_pos + overlap.dat.mer.k_len] = 0;
        getFragRecordSequence(dbgf2)[overlap.dat.mer.b_pos + overlap.dat.mer.k_len] = 0;

        fprintf(stderr, ">"F_IID"\n%s\n", getFragRecordIID(dbgf1), getFragRecordSequence(dbgf1) + overlap.dat.mer.a_pos);
        fprintf(stderr, ">"F_IID"\n%s\n", getFragRecordIID(dbgf2), getFragRecordSequence(dbgf2) + overlap.dat.mer.b_pos);

        assert(0);
      }
#endif

      //  Now, skip ahead until we find the next pair.
      //
      u64bit  lastiid = hits[i].tseq;
      while ((i < hitsLen) && (hits[i].tseq == lastiid))
        i++;

    }  //  over all hits

    delete M;
  }

  fprintf(stderr, "Found "u64bitFMT" mer hits.\n", merfound);
  fprintf(stderr, "Found "u64bitFMT" overlaps.\n", ovlfound);

  AS_OVS_closeBinaryOverlapFile(binout);

  closeFragStream(frs);
  closeGateKeeperStore(gkp);

  delete PS;
  delete MS;
  delete CS;

  delete gkpseq;
}
