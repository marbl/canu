
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-23 to 2015-JUL-26
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "sweatShop.H"
//#include <pthread.h>

#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "overlapReadCache.H"

#include "NDalign.H"
#include "analyzeAlignment.H"

#include "AS_UTL_reverseComplement.H"

#include "timeAndSize.H" //  getTime();



class consensusGlobalData {
public:
  consensusGlobalData(double  maxErate_,
                      uint32  bgnID_,
                      int32   endID_,
                      char   *gkpName,
                      char   *ovlName,
                      char   *tigName,  uint32  tigVers,
                      char   *cnsName,
                      char   *fastqName,
                      uint64 memLimit_) {

    //  Inputs

    gkpStore  = gkStore::gkStore_open(gkpName);

    readCache = new overlapReadCache(gkpStore, memLimit);

    ovlStore  = (ovlName) ? new ovStore(ovlName, gkpStore) : NULL;
    tigStore  = (tigName) ? new tgStore(tigName, tigVers)  : NULL;

    if (ovlStore)
      fprintf(stderr, "consensusGlobalData()--  opened ovlStore '%s'\n", ovlName);
    if (tigStore)
      fprintf(stderr, "consensusGlobalData()--  opened tigStore '%s'\n", tigName);

    //  Parameters

    maxErate   = maxErate_;
    memLimit   = memLimit_;

    if (bgnID_ == 0)                                bgnID_ = 1;
    if (endID_  > gkpStore->gkStore_getNumReads())  endID_ = gkpStore->gkStore_getNumReads() + 1;

    bgnID = bgnID_;
    curID = bgnID_;
    endID = endID_;

    //  Outputs

    cnsFile   = NULL;
    fastqFile = NULL;

    if (cnsName) {
      errno = 0;
      cnsFile = fopen(cnsName, "w");
      if (errno)
        fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", cnsName, strerror(errno)), exit(1);
    }

    if (fastqName) {
      errno = 0;
      fastqFile = fopen(fastqName, "w");
      if (errno)
        fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", fastqName, strerror(errno)), exit(1);
    }

    //  State for loading overlaps

    //  State for loading tigs

  };

  ~consensusGlobalData() {
    gkpStore->gkStore_close();

    delete readCache;
    delete ovlStore;
    delete tigStore;

    if (cnsFile)
      fclose(cnsFile);

    if (fastqFile)
      fclose(fastqFile);
  };

  //  Parameters

  double             maxErate;
  uint64             memLimit;

  uint32             bgnID;
  uint32             curID;  //  Currently loading id
  uint32             endID;

  //  Inputs

  gkStore           *gkpStore;

  overlapReadCache  *readCache;

  ovStore           *ovlStore;
  tgStore           *tigStore;

  //  State for loading

  gkReadData         readData;

  //  State for loading overlaps

  //  State for loading tigs

  //  Outputs

  FILE              *cnsFile;
  FILE              *fastqFile;
};



class consensusThreadData {
public:
  consensusThreadData(consensusGlobalData *g, uint32 tid) {
    threadID = tid;

    nPassed  = 0;
    nFailed  = 0;

    align    = new NDalign(pedGlobal, g->maxErate, 15);  //  true = partial aligns, maxErate, seedSize
    analyze  = new analyzeAlignment();
  };
  ~consensusThreadData() {
    delete align;
    delete analyze;
  };

  uint32                  threadID;

  uint64                  nPassed;
  uint64                  nFailed;

  char                    bRev[AS_MAX_READLEN];

  NDalign                *align;
  analyzeAlignment       *analyze;
};




//  The overlap compute needs both strings in the correct orientation.
//  This loaded just loads the tig/overlaps, and converts to a common format.
//  It ensures that the overlapReadCache is loaded.

//  aID aBgn aEnd
//  bID bBgn bEnd bFlip

class consensusComputation {
public:
  consensusComputation(tgTig  *tig) {
    _tig    = tig;
    _corLen = 0;
    _corSeq = NULL;
    _corQlt = NULL;
  };

  ~consensusComputation() {
    delete [] _corSeq;
    delete [] _corQlt;
  };


public:
  tgTig      *_tig;          //  Input

  //ovOverlap  *_overlaps;     //  Input, we convert this to a tig...
  //uint32      _overlapsLen;

  uint32      _corLen;
  char       *_corSeq;       //  Output sequence
  char       *_corQlt;
};




//  Complicated.  Load the overlaps, filter, and convert to a tig.
consensusComputation *
consensusReaderOverlaps(consensusGlobalData *UNUSED(g)) {
  return(NULL);
}


//  Simple, just load the tig and call it a day.
consensusComputation *
consensusReaderTigs(consensusGlobalData *g) {
  tgTig                 *t = NULL;
  consensusComputation  *s = NULL;

  while ((t == NULL) && (g->curID < g->endID))
    t = g->tigStore->loadTig(g->curID++);

  if (t)
    s = new consensusComputation(t);

  return(s);
}


void *
consensusReader(void *G) {
  consensusGlobalData    *g = (consensusGlobalData  *)G;
  consensusComputation   *s = NULL;

  if ((g->curID < g->endID) && (g->ovlStore))
    s = consensusReaderOverlaps(g);

  if ((g->curID < g->endID) && (g->tigStore))
    s = consensusReaderTigs(g);

  if (s)
    g->readCache->loadReads(s->_tig);

  return(s);
}


void
consensusWorker(void *G, void *T, void *S) {
  consensusGlobalData    *g = (consensusGlobalData  *)G;
  consensusThreadData    *t = (consensusThreadData  *)T;
  consensusComputation   *s = (consensusComputation *)S;

  uint32  rID = s->_tig->tigID();

  fprintf(stderr, "THREAD %u working on tig %u\n", t->threadID, rID);

  t->analyze->reset(rID,
                    g->readCache->getRead(rID),
                    g->readCache->getLength(rID));

  for (uint32 oo=0; oo<s->_tig->numberOfChildren(); oo++) {
    if (s->_tig->getChild(oo)->isRead() == false)
      continue;

    tgPosition  *pos = s->_tig->getChild(oo);

    //
    //  This colosely follows overlapPair
    //

    //  Load A.

    uint32  aID  = s->_tig->tigID();
    char   *aStr = g->readCache->getRead  (aID);
    uint32  aLen = g->readCache->getLength(aID);

    int32   aLo = pos->min() - 100;    if (aLo < 0)  aLo = 0;
    int32   aHi = pos->max() + 100;

    assert(aID == rID);
    assert(aLo < aHi);

    //  Load B.  If reversed, we need to reverse the coordinates to meet the overlap spec.

    uint32  bID  = pos->ident();
    char   *bStr = g->readCache->getRead  (bID);
    uint32  bLen = g->readCache->getLength(bID);

    int32   bLo = (pos->isReverse() == false) ? (       pos->askip()) : (bLen - pos->askip());
    int32   bHi = (pos->isReverse() == false) ? (bLen - pos->bskip()) : (       pos->bskip());

    //  Compute the overlap

#if 0
      fprintf(stderr, "aligned  %6u %s %6u -- %5u-%5u %5u-%5u\n",
              aID,
              pos->isReverse() ? "<--" : "-->",
              bID,
              aLo,  aHi,  bLo,  bHi);
#endif

    t->align->initialize(aID, aStr, aLen, aLo, aHi,
                         bID, bStr, bLen, bLo, bHi, pos->isReverse());

    if (t->align->findMinMaxDiagonal(40) == false)
      continue;

    if (t->align->findSeeds(false) == false)
      continue;

    t->align->findHits();
    t->align->chainHits();

    if (t->align->processHits() == true) {
      t->nPassed++;

      int32  aLoO = t->align->abgn();
      int32  aHiO = t->align->aend();

      int32  bLoO = t->align->bbgn();
      int32  bHiO = t->align->bend();

      if (bLoO > bHiO) {
        bLoO = t->align->bend();
        bHiO = t->align->bbgn();
      }

#if 0
      fprintf(stderr, "aligned  %6u %s %6u -- %5u-%5u %5u-%5u -- %5u-%5u %5u-%5u -- %6.2f\n",
              aID,
              pos->isReverse() ? "<--" : "-->",
              bID,
              aLo,  aHi,  bLo,  bHi,
              aLoO, aHiO, bLoO, bHiO,
              100.0 * (aHiO - aLoO) / (aHi - aLo));
      //t->align->display();
#endif

      //  This wants pointers to the start of the strings that align, and the offset to the full read.
      t->analyze->analyze(t->align->astr() + aLoO, aHiO - aLoO, aLoO,
                          t->align->bstr() + bLoO, bHiO - bLoO,
                          t->align->deltaLen(),
                          t->align->delta());
    } else {
      t->nFailed++;

#if 1
      fprintf(stderr, "FAILED   %6u %s %6u -- %5u-%5u %5u-%5u -- %5u-%5u %5u-%5u -- %6.2f\n",
              aID,
              pos->isReverse() ? "<--" : "-->",
              bID,
              aLo,  aHi,  bLo,  bHi,
              0, 0, 0, 0,
              0.0);
      //t->align->display();
#endif
    }
  }

  //fprintf(stderr, "THREAD %2u finished with tig %5u -- passed %12u -- failed %12u\n", t->threadID, s->_tig->tigID(), t->nPassed, t->nFailed);

  t->analyze->generateCorrections();
  t->analyze->generateCorrectedRead();

  s->_corLen = t->analyze->_corSeqLen;
  s->_corSeq = new char [s->_corLen + 1];
  s->_corQlt = new char [s->_corLen + 1];

  memcpy(s->_corSeq, t->analyze->_corSeq, sizeof(char) * (s->_corLen + 1));
  memcpy(s->_corQlt, t->analyze->_corSeq, sizeof(char) * (s->_corLen + 1));

  for (uint32 ii=0; ii<s->_corLen; ii++)
    s->_corQlt[ii] = '!' + 10;
}


void
consensusWriter(void *G, void *S) {
  consensusGlobalData    *g = (consensusGlobalData  *)G;
  consensusComputation   *s = (consensusComputation *)S;

  if (g->cnsFile) {
  }

  if (g->fastqFile) {
    fprintf(g->fastqFile, "@read%08u\n", s->_tig->tigID());
    fprintf(g->fastqFile, "%s\n", s->_corSeq);
    fprintf(g->fastqFile, "+\n");
    fprintf(g->fastqFile, "%s\n", s->_corQlt);
  }

  delete s;
}







int
main(int argc, char **argv) {
  char    *gkpName         = NULL;
  char    *ovlName         = NULL;
  char    *tigName         = NULL;
  uint32   tigVers         = 1;

  char    *cnsName         = NULL;
  char    *fastqName       = NULL;

  uint32   bgnID           = 1;
  uint32   endID           = UINT32_MAX;

  uint32   numThreads      = 1;

  double   maxErate        = 0.02;
  uint64   memLimit        = 4;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      cnsName = argv[++arg];
    } else if (strcmp(argv[arg], "-f") == 0) {
      fastqName = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-t") == 0) {
      numThreads = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-erate") == 0) {
      maxErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memLimit = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if ((ovlName == NULL) && (tigName == NULL))
    err++;
  if ((ovlName != NULL) && (tigName != NULL))
    err++;
  if ((cnsName == NULL) && (fastqName == NULL))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -G gkpStore     Mandatory, path to gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Inputs can come from either an overlap or a tig store.\n");
    fprintf(stderr, "  -O ovlStore     \n");
    fprintf(stderr, "  -T tigStore tigVers      \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If from an ovlStore, the range of reads processed can be restricted.\n");
    fprintf(stderr, "  -b bgnID        \n");
    fprintf(stderr, "  -e endID        \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs will be written as the full multialignment and the final consensus sequence\n");
    fprintf(stderr, "  -c output.cns   \n");
    fprintf(stderr, "  -f output.fastq \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -erate e        Overlaps are computed at 'e' fraction error; must be larger than the original erate\n");
    fprintf(stderr, "  -memory m       Use up to 'm' GB of memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t n            Use up to 'n' cores\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no gatekeeper (-G) supplied.\n");
    if ((ovlName == NULL) && (tigName == NULL))
      fprintf(stderr, "ERROR: no inputs (-O or -T) supplied.\n");
    if ((ovlName != NULL) && (tigName != NULL))
      fprintf(stderr, "ERROR: only one input (-O or -T) may be supplied.\n");
    if ((cnsName == NULL) && (fastqName == NULL))
      fprintf(stderr, "ERROR: no outputs (-c or -f) supplied.\n");

    exit(1);
  }

  consensusGlobalData  *g = new consensusGlobalData(maxErate,
                                                    bgnID,
                                                    endID,
                                                    gkpName,
                                                    ovlName,
                                                    tigName, tigVers,
                                                    cnsName,
                                                    fastqName,
                                                    memLimit);

#if 1

  consensusThreadData  *t = new consensusThreadData(g, 0);

  while (1) {
    consensusComputation *c = (consensusComputation *)consensusReader(g);

    if (c == NULL)
      break;

    consensusWorker(g, t, c);
    consensusWriter(g, c);
  }

  delete t;

#else

  consensusThreadData **td = new consensusThreadData * [numThreads];
  sweatShop            *ss = new sweatShop(consensusReader, consensusWorker, consensusWriter);

  ss->setLoaderQueueSize(16384);
  ss->setWriterQueueSize(1024);

  ss->setNumberOfWorkers(numThreads);

  for (uint32 w=0; w<numThreads; w++)
    ss->setThreadData(w, td[w] = new consensusThreadData(g, w));  //  these leak

  ss->run(g, true);

  delete ss;

  for (uint32 w=0; w<numThreads; w++)
    delete td[w];

#endif

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}
