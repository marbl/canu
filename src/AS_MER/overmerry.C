
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
#include "sweatShop.H"


//  using external mer counts - we'll need to extend the posDB to have
//  more stuff in it, or construct another mer lookup hash.
//
//  parallel searching
//
//  subset of gkpstore
//
//  extend merstream with spaced, compressed, skips, transitions




//  Instead of using The internal overlap, which has enough extra
//  stuff in it that we cannot store a sequence iid for the table
//  sequence, we need to make our own overlap structure.
//
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





class ovmGlobalData {
public:
  ovmGlobalData() {
    gkpPath     = 0L;
    merSize     = 23;
    compression = 1;
    maxCount    = 40;
    numThreads  = 4;

    qGK  = 0L;
    qFS  = 0L;
    qBeg = 1;
    qEnd = 0;

    tGK  = 0L;
    tKB  = 0L;
    tSS  = 0L;
    tMS  = 0L;
    tPS  = 0L;
    tBeg = 1;
    tEnd = 0;

    outputName  = 0L;
    outputFile = 0L;
  };

  ~ovmGlobalData() {

    //fprintf(stderr, "Found "u64bitFMT" mer hits.\n", merfound);
    //fprintf(stderr, "Found "u64bitFMT" overlaps.\n", ovlfound);

    AS_OVS_closeBinaryOverlapFile(outputFile);

    closeFragStream(qFS);
    closeGateKeeperStore(qGK);

    delete tPS;
    delete tMS;
    delete tSS;
    delete tKB;
    delete tGK;
  };

  void  build(void) {

    //
    //  Open inputs for the reader
    //
    qGK = openGateKeeperStore(gkpPath, FALSE);

    //  Use that gkpStore to check and reset the desired ranges
    //
    if (qBeg < 1)   qBeg = 1;
    if (qEnd == 0)  qEnd = getNumGateKeeperFragments(qGK);
    if (qBeg >= qEnd) {
      fprintf(stderr, "ERROR: -qb="F_U32" and -qe="F_U32" are invalid ("F_U32" frags in the store)\n",
              qBeg, qEnd, getNumGateKeeperFragments(qGK));
      exit(1);
    }

    if (tBeg < 1)   tBeg = 1;
    if (tEnd == 0)  tEnd = getNumGateKeeperFragments(qGK);
    if (tBeg >= tEnd) {
      fprintf(stderr, "ERROR: -tb="F_U32" and -te="F_U32" are invalid ("F_U32" frags in the store)\n",
              tBeg, tEnd, getNumGateKeeperFragments(qGK));
      exit(1);
    }

    //  Use that gkpStore to quickly build a list of the clear ranges
    //  and full sequence length for reads in the table.
    //
    {
      fragRecord   fr;
      FragStream  *fs = openFragStream(qGK, FRAG_S_INF);

      resetFragStream(fs, tBeg, tEnd);

      table_clrBeg       = new uint32 [tEnd - tBeg];
      table_untrimLength = new uint32 [tEnd - tBeg];

      while (nextFragStream(fs, &fr)) {
        table_clrBeg      [getFragRecordIID(&fr) - tBeg] = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBT);
        table_untrimLength[getFragRecordIID(&fr) - tBeg] = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_UNTRIM);
      }

      closeFragStream(fs);
    }

    //  Open another fragStream for the query fragments.
    //
    qFS = openFragStream(qGK, FRAG_S_SEQ);
    resetFragStream(qFS, qBeg, qEnd);


    //
    //  Build state for the workers
    //

    //  Open the gatekeeper store as a kmer seqFile.  We need to
    //  dynamic_cast this back to our gkpStoreSequence, so we can access
    //  methods defined only on that object.
    //
    {
      char     gkpName[FILENAME_MAX + 64] = {0};
      sprintf(gkpName, "%s:%u-%u:obt", gkpPath, tBeg, tEnd);

      tGK = dynamic_cast<gkpStoreSequence*>(openSeqFile(gkpName));
      if (tGK == 0L) {
        fprintf(stderr, "%s: invalid input file '%s' (not a GateKeeperStore?).\n", gkpName);
        exit(1);
      }

      tKB = new kMerBuilder(merSize, compression, 0L);

      tSS = new seqStream(tGK, true);
      tSS->setSeparator('.', 1);

      tMS = new merStream(tKB, tSS);
      tPS = new positionDB(tMS, merSize, 0, 0L, 0L, maxCount, true);  //  This interface is in kmer r1544
    }

    //
    //  Open the output file.
    //

    outputFile = AS_OVS_createBinaryOverlapFile(outputName, FALSE);
  };

  uint32    getClrBeg(CDS_IID_t iid) {
    assert(tBeg <= iid);
    return(table_clrBeg[iid - tBeg]);
  };

  uint32    getUntrimLength(CDS_IID_t iid) {
    assert(tBeg <= iid);
    return(table_untrimLength[iid - tBeg]);
  };


  //  Command line parameters
  //
  char    *gkpPath;
  uint32   merSize;
  uint32   compression;
  uint32   maxCount;
  uint32   numThreads;

  //  for the READER only
  //
  GateKeeperStore   *qGK;
  FragStream        *qFS;
  uint32             qBeg;
  uint32             qEnd;

  //  for the WORKERS.
  //
  gkpStoreSequence  *tGK;
  kMerBuilder       *tKB;
  seqStream         *tSS;  //  needs to be public so we can offset coords
  merStream         *tMS;
  positionDB        *tPS;  //  needs to be public!  (this is the main tabile)
  uint32             tBeg;
  uint32             tEnd;

  //  for the WORKERS - we need to know the clrBeg and full read
  //  length for stuff in the table, which we can't get from the
  //  READER gkp.
  //
  uint32            *table_clrBeg;
  uint32            *table_untrimLength;

  //  for the WRITER only.
  //
  char              *outputName;
  BinaryOverlapFile *outputFile;
};



class ovmThreadData {
public:
  ovmThreadData(ovmGlobalData *g) {
    qKB      = new kMerBuilder(g->merSize, g->compression, 0L);

    posnF    = 0L;
    posnFMax = 0;
    posnFLen = 0;

    posnR    = 0L;
    posnRMax = 0;
    posnRLen = 0;

    hitsLen = 0;
    hitsMax = 0;
    hits    = 0L;

    merfound = 0;
    ovlfound = 0;
  };

  ~ovmThreadData() {
    delete [] posnF;
    delete [] posnR;
    delete [] hits;
  };

  void
  addHit(seqStream   *SS,
         CDS_IID_t    iid,
         u64bit       qpos,
         u64bit       pos,
         u64bit       cnt,
         u64bit       pal,
         u64bit       fwd) {

    uint32  seq = SS->sequenceNumberOfPosition(pos);

    pos -= SS->startOf(seq);
    seq  = SS->IIDOf(seq);

    if (iid == seq)
      return;

    if (iid < seq)
      return;

    if (hitsLen >= hitsMax) {
      if (hitsMax == 0) {
        hitsMax = 1048576;  //  tiny, 8MB per thread
        hits    = new kmerhit [hitsMax];
      } else {
        hitsMax *= 2;
        kmerhit *h = new kmerhit [hitsMax];
        memcpy(h, hits, sizeof(kmerhit) * hitsLen);
        delete [] hits;
        hits = h;
      }
    }

    hits[hitsLen].tseq = seq;
    hits[hitsLen].tpos = pos;
    hits[hitsLen].qpos = qpos;
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
    merfound++;
  };


  kMerBuilder  *qKB;

  u64bit        posnFLen;
  u64bit        posnFMax;
  u64bit       *posnF;

  u64bit        posnRLen;
  u64bit        posnRMax;
  u64bit       *posnR;

  u32bit        hitsLen;
  u32bit        hitsMax;
  kmerhit      *hits;

  u64bit        merfound;
  u64bit        ovlfound;
};



class ovmComputation {
public:
  ovmComputation(fragRecord *fr) {
    beg = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBT);
    end = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBT);
    tln = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_UNTRIM);

    iid = getFragRecordIID(fr);
    uid = getFragRecordUID(fr);

    memset(seq, 0, AS_FRAG_MAX_LEN);
    strcpy(seq, getFragRecordSequence(fr));

    ovsLen = 0;
    ovsMax = 1024;
    ovs    = new OVSoverlap [ovsMax];
  };

  ~ovmComputation() {
    delete [] ovs;
  };

  void        addOverlap(OVSoverlap *overlap) {

#if 0
    //  some ascii output
    fprintf(stderr, "%d\t%d\t%c\t%d\t%d\n",
            (int)overlap->a_iid,
            (int)overlap->dat.mer.a_pos,
            overlap->dat.mer.palindrome ? 'p' : (overlap->dat.mer.fwd ? 'f' : 'r'),
            (int)overlap->b_iid,
            (int)overlap->dat.mer.b_pos);
#endif

    if (ovsLen >= ovsMax) {
      ovsMax *= 2;
      OVSoverlap *o = new OVSoverlap [ovsMax];
      memcpy(o, ovs, sizeof(OVSoverlap) * ovsLen);
      delete [] ovs;
      ovs = o;
    }

    ovs[ovsLen++] = *overlap;
  };

  void        writeOverlaps(BinaryOverlapFile *outputFile) {
    for (uint32 i=0; i<ovsLen; i++)
      AS_OVS_writeOverlap(outputFile, ovs + i);
  };

  uint32      beg;
  uint32      end;
  uint32      tln;

  CDS_IID_t   iid;
  CDS_UID_t   uid;

  uint32      ovsLen;  //  Overlap Storage, waiting for output
  uint32      ovsMax;
  OVSoverlap *ovs;

  char        seq[AS_FRAG_MAX_LEN];
};




void*
ovmReader(void *G) {
  ovmGlobalData  *g = (ovmGlobalData *)G;

  static fragRecord  fr;  //  static only for performance

 again:
  if (nextFragStream(g->qFS, &fr) == 0)
    return(0L);

  if (getFragRecordIsDeleted(&fr))
    goto again;

  return(new ovmComputation(&fr));
}




void
ovmWorker(void *G, void *T, void *S) {
  ovmGlobalData    *g = (ovmGlobalData  *)G;
  ovmThreadData    *t = (ovmThreadData  *)T;
  ovmComputation   *s = (ovmComputation *)S;

  OVSoverlap        overlap;

#if 0
  if (s->iid == 76)
    fprintf(stderr, "IID 76\n");
#endif

  t->hitsLen = 0;

  //  We can "simulate" a canonical mercount by getting the count
  //  for both forward and reverse.

  merStream *sMSTR  = new merStream(t->qKB, s->seq, s->beg, s->end - s->beg);
  uint32    *sSPAN  = new uint32 [s->end - s->beg];

  while (sMSTR->nextMer()) {
    u64bit  qpos = sMSTR->thePositionInSequence();

    sSPAN[qpos] = sMSTR->theFMer().getMerSpan();
    assert(qpos <= s->end - s->beg);

    if (sMSTR->theFMer() == sMSTR->theRMer()) {
      g->tPS->get(sMSTR->theFMer(), t->posnF, t->posnFMax, t->posnFLen);
      if (t->posnFLen > 1)
        for (u32bit i=0; i<t->posnFLen; i++)
          t->addHit(g->tSS, s->iid, qpos, t->posnF[i], t->posnFLen, 1, 0);
    } else {
      g->tPS->get(sMSTR->theFMer(), t->posnF, t->posnFMax, t->posnFLen);
      g->tPS->get(sMSTR->theRMer(), t->posnR, t->posnRMax, t->posnRLen);

      u64bit totalLen = t->posnFLen + t->posnRLen - 1;  //  the canonical mer count

      if (t->posnFLen > 1)
        for (u32bit i=0; i<t->posnFLen; i++)
          t->addHit(g->tSS, s->iid, qpos, t->posnF[i], totalLen, 0, 1);
      if (t->posnRLen > 1)
        for (u32bit i=0; i<t->posnRLen; i++)
          t->addHit(g->tSS, s->iid, qpos, t->posnR[i], totalLen, 0, 0);
    }
  }

  delete sMSTR;

  if (t->hitsLen == 0)
    return;


  //  We have all the hits for this frag.  Sort them by sequence
  //  (the other sequence), then pick out the one with the least
  //  count for each sequence.

  qsort(t->hits, t->hitsLen, sizeof(kmerhit), kmerhitcompare);


#if 0
  //  Debug, I guess.  Generates lots of output, since frags with
  //  big identical overlaps will have lots and lots of mers in
  //  common.
  //
  for (u32bit i=0; i<t->hitsLen; i++) {
    if (i != t->hitsLen) {
      fprintf(stderr, u32bitFMT"\t"u64bitFMT"\t"u32bitFMT"\t"u64bitFMT"\t%c\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\tTAG\n",
              t->hits[i].tseq, t->hits[i].tpos,
              s->iid,  t->hits[i].qpos,
              t->hits[i].pal ? 'p' : (t->hits[i].fwd ? 'f' : 'r'),
              0,
              t->hits[i].cnt,
              merSize);
    }
  }
#endif



  for (u32bit i=0; i<t->hitsLen; i++) {

#if 0
    if ((t->hits[i].tseq == 48) && (s->iid == 76))
      fprintf(stderr, "Found it.\n");
#endif

    //  By the definition of our sort, the least common mer is the
    //  first hit in the list for each pair of sequences.
    //
    t->ovlfound++;

    //  Adjust coords to be relative to whole read -- the table is
    //  built using only sequence in the OBT clear.  The query is
    //  built starting at the begin of the OBT clear.  Same effect for
    //  both just a different mechanism.
    //
    t->hits[i].tpos += g->getClrBeg(t->hits[i].tseq);
    t->hits[i].qpos += s->beg;

    //  Reverse if needed -- we need to remember, from when we were
    //  grabbing mers, the length of the uncompressed mer -- that's
    //  the sSPAN; if we're not using compressed seeds, sSPAN ==
    //  g->merSize.  The [t->hits[i].qpos - s->beg] array index is
    //  simply the position in the trimmed read.
    //
    if (t->hits[i].fwd == false)
      t->hits[i].qpos = s->tln - t->hits[i].qpos - sSPAN[t->hits[i].qpos - s->beg];

    //  Save off the A vs B overlap
    //
    overlap.a_iid                      = t->hits[i].tseq;
    overlap.b_iid                      = s->iid;
    overlap.dat.mer.datpad             = 0;
    overlap.dat.mer.compression_length = g->compression;
    overlap.dat.mer.fwd                = t->hits[i].fwd;
    overlap.dat.mer.palindrome         = t->hits[i].pal;
    overlap.dat.mer.a_pos              = t->hits[i].tpos;
    overlap.dat.mer.b_pos              = t->hits[i].qpos;
    overlap.dat.mer.k_count            = t->hits[i].cnt;
    overlap.dat.mer.k_len              = g->merSize;
    overlap.dat.mer.type               = AS_OVS_TYPE_MER;
    s->addOverlap(&overlap);

    //  Save off the B vs A overlap
    //
    overlap.a_iid = s->iid;
    overlap.b_iid = t->hits[i].tseq;

    if (overlap.dat.mer.fwd) {
      overlap.dat.mer.a_pos = t->hits[i].qpos;
      overlap.dat.mer.b_pos = t->hits[i].tpos;
    } else {
      uint32 othlen = g->getUntrimLength(t->hits[i].tseq);

      //  The -1 is to back up to the last base in the mer.

      overlap.dat.mer.a_pos = s->tln - t->hits[i].qpos - 1;
      overlap.dat.mer.b_pos = othlen - t->hits[i].tpos - 1;

#if ADJUST_LEFT_END
      //  This only works for non-compressed seeds.
      overlap.dat.mer.a_pos -= g->merSize;  //  sSPAN[t->hits[i].qpos];
      overlap.dat.mer.b_pos -= g->merSize;  //  we don't have the span of the table-based mer
#endif

    }
    s->addOverlap(&overlap);

    //  Now, skip ahead until we find the next pair.
    //
    u64bit  lastiid = t->hits[i].tseq;
    while ((i < t->hitsLen) && (t->hits[i].tseq == lastiid))
      i++;
  }  //  over all hits
}



void
ovmWriter(void *G, void *S) {
  ovmGlobalData    *g = (ovmGlobalData  *)G;
  ovmComputation   *s = (ovmComputation *)S;

  s->writeOverlaps(g->outputFile);

  delete s;
}





int
main(int argc, char **argv) {
  ovmGlobalData  *g = new ovmGlobalData;

  assert(sizeof(kmerhit) == 8);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      g->gkpPath = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      g->merSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      g->compression = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      g->outputName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-tb") == 0) {
      g->tBeg = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-te") == 0) {
      g->tEnd = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-qb") == 0) {
      g->qBeg = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-qe") == 0) {
      g->qEnd = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((g->gkpPath == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-g gkpStore] [-m merSize] [-c compression] [-o outputName]\n", argv[0]);
    exit(1);
  }

  seqFactory::instance()->registerFile(new gkpStoreSequence());

  g->build();

  sweatShop *ss = new sweatShop(ovmReader, ovmWorker, ovmWriter);

  ss->loaderQueueSize(10240);
  ss->writerQueueSize(10240);

  ss->numberOfWorkers(g->numThreads);

  for (u32bit w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new ovmThreadData(g));

  ss->run(g, true);  //  true == verbose

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");
}
