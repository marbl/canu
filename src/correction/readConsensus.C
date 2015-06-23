const char *mainid = "$Id:  $";

#include "AS_global.H"

#include "sweatShop.H"
//#include <pthread.h>

#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "overlapReadCache.H"

#include "overlapAlign.H"
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

    //  Parameters

    maxErate   = maxErate_;
    memLimit   = memLimit_;

    bgnID = bgnID_;
    curID = bgnID_;
    endID = endID_;

    //  Inputs

    gkpStore  = new gkStore(gkpName);

    readCache = new overlapReadCache(gkpStore, memLimit);

    ovlStore  = (ovlName) ? new ovStore(ovlName, gkpStore) : NULL;
    tigStore  = (tigName) ? new tgStore(tigName, tigVers)  : NULL;

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
    delete gkpStore;
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
    align    = new overlapAlign(true, g->maxErate, 15);  //  partial aligns, maxErate, seedSize
    analyze  = new analyzeAlignment();
  };
  ~consensusThreadData() {
    delete align;
    delete analyze;
  };

  uint32                 threadID;

  char                    bRev[AS_MAX_READLEN];

  overlapAlign           *align;
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
    _corCns = NULL;
    _corQlt = NULL;
    _corLen = NULL;
  };

  ~consensusComputation() {
    delete [] _corCns;
    delete [] _corQlt;
    delete [] _corLen;
  };
  

public:
  tgTig      *_tig;          //  Input

  ovOverlap  *_overlaps;     //  Input, we convert this to a tig...
  uint32      _overlapsLen;

  char       *_corCns;       //  Output sequence
  char       *_corQlt;
  char       *_corLen;
};




//  Complicated.  Load the overlaps, filter, and convert to a tig.
consensusComputation *
consensusReaderOverlaps(consensusGlobalData *g) {
}


//  Simple, just load the tig and call it a day.
consensusComputation *
consensusReaderTigs(consensusGlobalData *g) {
  return(new consensusComputation(g->tigStore->loadTig(g->curID++)));
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

  fprintf(stderr, "WORKING on tig %u\n", s->_tig->tigID());

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

    int32   aLo = pos->min() - 100;
    int32   aHi = pos->max() + 100;

    assert(aLo < aHi);

    //  Load B.

    uint32  bID  = pos->ident();
    char   *bStr = g->readCache->getRead  (bID);
    uint32  bLen = g->readCache->getLength(bID);

    int32   bLo =        pos->askip();
    int32   bHi = bLen - pos->bskip();

    //  Make B the correct orientation, and adjust coordinates.

#if 0
    if (pos->isReverse()) {
      memcpy(t->bRev, bStr, sizeof(char) * (bLen + 1));

      reverseComplementSequence(t->bRev, bLen);

      bStr = t->bRev;

      //bLo =
      //bHi =
    }

    assert(bLo < bHi);
#endif

    //  Compute the overlap

    t->align->initialize(aID, aStr, aLen, aLo, aHi,
                         bID, bStr, bLen, bLo, bHi, pos->isReverse());

    if (t->align->findMinMaxDiagonal(40) == false) {
      //fprintf(stderr, "A %6u %5d-%5d ->   B %6u %5d-%5d %s ALIGN LENGTH TOO SHORT.\n",
      //        aID, ovl->a_bgn(), ovl->a_end(),
      //        bID, ovl->b_bgn(), ovl->b_end(),
      //        ovl->flipped() ? "<-" : "->");
      continue;
    }

    if (t->align->findSeeds(false) == false) {
      //fprintf(stderr, "A %6u %5d-%5d ->   B %6u %5d-%5d %s NO SEEDS.\n",
      //        aID, ovl->a_bgn(), ovl->a_end(),
      //        bID, ovl->b_bgn(), ovl->b_end(),
      //        ovl->flipped() ? "<-" : "->");
      continue;
    }

    t->align->findHits();
    t->align->chainHits();


#if 0
    if (t->align->processHits() == true) {
      //analyze->analyze(aStr, aLen, aOffset,
      //                 bStr, bLen,
      //                 deltaLen,
      //                 delta);
    } else {
    }
#endif
  }

  //analyze->generateCorrections(NULL);
  //analyze->generateCorrectedRead(NULL, NULL, NULL);

  //
  //  END of following
  //



}


void
consensusWriter(void *G, void *S) {
  consensusGlobalData    *g = (consensusGlobalData  *)G;
  consensusComputation   *s = (consensusComputation *)S;

  if (g->cnsFile) {
  }

  if (g->fastqFile) {
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

  uint32   bgnID           = 0;
  uint32   endID           = UINT32_MAX;

  uint32   numThreads      = 1;

  double   maxErate        = 0.12;
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
  if ((ovlName == NULL) || (tigName == NULL))
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
    if ((ovlName == NULL) || (tigName == NULL))
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

  sweatShop  *ss = new sweatShop(consensusReader, consensusWorker, consensusWriter);

  ss->setLoaderQueueSize(16384);
  ss->setWriterQueueSize(1024);

  ss->setNumberOfWorkers(numThreads);

  for (uint32 w=0; w<numThreads; w++)
    ss->setThreadData(w, new consensusThreadData(g, w));  //  these leak

  ss->run(g, false);

  delete ss;
  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}

 

