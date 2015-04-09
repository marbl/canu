const char *mainid = "$Id:  $";

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "stashContains.H"


//  Generate a layout for the read in ovl[0].a_iid, using most or all of the overlaps
//  in ovl.

tgTig *
generateLayout(gkStore    *gkpStore,
               uint32      minLen,
               double      maxErate,
               double      maxCoverage,
               ovsOverlap *ovl,
               uint32      ovlLen) {
  //fprintf(stderr, "Generate layout for read "F_U32" with "F_U32" overlaps.\n",
  //        ovl[0].a_iid, ovlLen);

  tgTig  *layout = new tgTig;

  layout->_tigID           = ovl[0].a_iid;
  layout->_coverageStat    = 1.0;  //  Default to just barely unique
  layout->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  layout->_suggestRepeat   = false;
  layout->_suggestUnique   = false;
  layout->_suggestCircular = false;
  layout->_suggestHaploid  = false;

  gkRead  *read = gkpStore->gkStore_getRead(ovl[0].a_iid);

  layout->_layoutLen       = read->gkRead_sequenceLength();

  resizeArray(layout->_children, layout->_childrenLen, layout->_childrenMax, ovlLen, resizeArray_doNothing);

  for (uint32 oo=0; oo<ovlLen; oo++) {
    if (ovl[oo].erate() > maxErate) {
      //fprintf(stderr, "  filter read %9u at position %6u,%6u erate %.3f - low quality (threshold %.2f)\n",
      //        ovl[oo].b_iid, ovl[oo].a_bgn(gkpStore), ovl[oo].a_end(gkpStore), ovl[oo].erate(), maxErate);
      continue;
    }

    if (ovl[oo].a_end(gkpStore) - ovl[oo].a_bgn(gkpStore) < minLen) {
      //fprintf(stderr, "  filter read %9u at position %6u,%6u erate %.3f - too short (threshold %u)\n",
      //        ovl[oo].b_iid, ovl[oo].a_bgn(gkpStore), ovl[oo].a_end(gkpStore), ovl[oo].erate(), minLen);
      continue;
    }

    tgPosition   *pos = layout->addChild();

    pos->_objID       = ovl[oo].b_iid;
    pos->_isRead      = true;
    pos->_isUnitig    = false;
    pos->_isContig    = false;

    pos->_askip       = ovl[oo].dat.ovl.bhg5;
    pos->_bskip       = ovl[oo].dat.ovl.bhg3;

    //  Set the read.  Parent is always the read we're building for, hangs and position come from
    //  the overlap.  Easy as pie!

    //  Well, except the hangs are abused here to indicate the amount of the read that
    //  doesn't align.

    if (ovl[oo].flipped() == false) {
      pos->set(ovl[oo].a_iid,
               ovl[oo].a_hang(),
               ovl[oo].b_hang(),
               ovl[oo].a_bgn(gkpStore), ovl[oo].a_end(gkpStore));

    } else {
      pos->set(ovl[oo].a_iid,
               ovl[oo].a_hang(),
               ovl[oo].b_hang(),
               ovl[oo].a_end(gkpStore), ovl[oo].a_bgn(gkpStore));
    }
  }

  //  Use utgcns's stashContains to get rid of extra coverage; we don't care about it, and
  //  just delete it immediately.

  savedChildren *sc = stashContains(layout, maxCoverage, false);

  if (sc) {
    delete sc->children;
    delete sc;
  }

  //  stashContains also sorts by position, so we're done.

#if 0
  for (uint32 ii=0; ii<layout->numberOfChildren(); ii++) {
    tgPosition *pos = layout->getChild(ii);

    fprintf(stderr, "  read %9u at position %6u,%6u hangs %6d %6d %c unAl %5d %5d\n",
            pos->_objID,
            pos->_min, pos->_max,
            pos->_ahang, pos->_bhang,
            pos->isForward() ? 'F' : 'R',
            pos->_askip, pos->_bskip);
  }
#endif

  return(layout);
}






int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *ovlName   = 0L;
  char             *tigName   = 0L;

  uint32            errorRate = AS_OVS_encodeQuality(0.015);

  char             *outputPrefix  = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;

  argc = AS_configure(argc, argv);

  uint32            minEvidenceOverlap  = 40;
  uint32            minEvidenceCoverage = 1;

  uint32            iidMin   = 0;
  uint32            iidMax   = UINT32_MAX;

  uint32            minLen       = 0;
  double            maxErate     = 1.0;
  double            maxCoverage  = DBL_MAX;

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax  = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-L") == 0) {
      minLen  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      maxErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-C") == 0) {
      maxCoverage = atof(argv[++arg]);


    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (err) {
    exit(1);
  }

  gkStore  *gkpStore = new gkStore(gkpName);
  ovStore  *ovlStore = new ovStore(ovlName);
  tgStore  *tigStore = new tgStore(tigName);

  ovlStore->setRange(iidMin, iidMax);

  uint32       ovlMax = 1024 * 1024;
  uint32       ovlLen = 0;
  ovsOverlap  *ovl    = new ovsOverlap [ovlMax];

  ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true);

  while (ovlLen > 0) {
    tgTig *tig = generateLayout(gkpStore, minLen, maxErate, maxCoverage, ovl, ovlLen);

    tigStore->insertTig(tig, false);

    ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true);

    delete tig;
  }

  delete tigStore;
  delete ovlStore;
  delete gkpStore;

  return(0);
}

