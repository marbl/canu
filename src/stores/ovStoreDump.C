
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
 *  This file is derived from:
 *
 *    src/AS_OVS/overlapStore.C
 *    src/AS_OVS/overlapStore.c
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-MAR-08 to 2013-AUG-01
 *      are Copyright 2007-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-DEC-16 to 2009-AUG-14
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gregory Sims from 2012-FEB-01 to 2012-MAR-14
 *      are Copyright 2012 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-22 to 2015-JUN-25
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_UTL_decodeRange.H"
#include "splitToWords.H"

#include "gkStore.H"
#include "ovStore.H"



class dumpParameters {
public:
  dumpParameters() {
    no5p               = false;
    no3p               = false;
    noContainer        = false;
    noContained        = false;
    noRedundant        = false;

    noBogartContained  = false;
    noBogartSuspicious = false;
    noBogartSingleton  = false;

    erateMin           = 0.0;
    erateMax           = 1.0;

    lengthMin          = 0;
    lengthMax          = UINT32_MAX;

    queryMin           = 0;
    queryMax           = UINT32_MAX;


    ovlKept            = 0;
    ovlFiltered        = 0;

    ovl5p              = 0;
    ovl3p              = 0;
    ovlContainer       = 0;
    ovlContained       = 0;
    ovlRedundant       = 0;

    ovlErateLo         = 0;
    ovlErateHi         = 0;

    ovlLengthLo        = 0;
    ovlLengthHi        = 0;
  };

  bool        parametersAreDefaults(void) {       //  Returns true if all the parameters
    return((no5p               == false) &&       //  are the default setting.
           (no3p               == false) &&       //
           (noContainer        == false) &&       //  This looks horrid, but it was easier
           (noContained        == false) &&       //  and simpler than tracking what
           (noRedundant        == false) &&       //  was set.  I'd need to make setter
           (noBogartContained  == false) &&       //  functions for everything, make the
           (noBogartSuspicious == false) &&       //  setter set a 'isChanged' flag, and
           (noBogartSingleton  == false) &&       //  then update command line parsing.
           (erateMin           == 0.0) &&         //
           (erateMax           == 1.0) &&         //  It took longer to write this comment
           (lengthMin          == 0) &&           //  than to write this function.
           (lengthMax          == UINT32_MAX) &&  //
           (queryMin           == 0) &&           //  So there.
           (queryMax           == UINT32_MAX));
  };

  bool        filterOverlap(ovOverlap *overlap) {
    double erate    = overlap->erate();
    uint32 length   = (lengthMax - lengthMin) / 2;   //  until we compute it
    int32  ahang    = overlap->a_hang();
    int32  bhang    = overlap->b_hang();
    bool   filtered = false;

    if ((no5p == true) && (ahang < 0) && (bhang < 0)) {
      ovl5p++;
      filtered = true;
    }

    if ((no3p == true) && (ahang > 0) && (bhang > 0)) {
      ovl3p++;
      filtered = true;
    }

    if ((noContainer) && (ahang <= 0) && (bhang >= 0)) {
      ovlContainer++;
      filtered = true;
    }

    if ((noContained) && (ahang >= 0) && (bhang <= 0)) {
      ovlContained++;
      filtered = true;
    }

    if ((noRedundant) && (overlap->a_iid >= overlap->b_iid)) {
      ovlRedundant++;
      filtered = true;
    }

    if ((overlap->b_iid < queryMin) ||
        (queryMax < overlap->b_iid)) {
      filtered = true;
    }

    if (erate < erateMin) {
      ovlErateLo++;
      filtered = true;
    }

    if (erate > erateMax) {
      ovlErateHi++;
      filtered = true;
    }

    if (length < lengthMin) {
      ovlLengthLo++;
      filtered = true;
    }

    if (length > lengthMax) {
      ovlLengthHi++;
      filtered = true;
    }

    return(filtered);
  };

  //  Filtering options.

  bool        no5p;
  bool        no3p;
  bool        noContainer;
  bool        noContained;
  bool        noRedundant;

  bool        noBogartContained;
  bool        noBogartSuspicious;
  bool        noBogartSingleton;

  uint32      queryMin;
  uint32      queryMax;

  double      erateMin;
  double      erateMax;

  uint32      lengthMin;
  uint32      lengthMax;

  //  Counts of what we filtered.

  uint64      ovlKept;
  uint64      ovlFiltered;

  uint64      ovl5p;
  uint64      ovl3p;
  uint64      ovlContainer;
  uint64      ovlContained;
  uint64      ovlRedundant;

  uint64      ovlErateHi;
  uint64      ovlErateLo;

  uint64      ovlLengthHi;
  uint64      ovlLengthLo;
};



int
main(int argc, char **argv) {
  char                 *gkpName     = NULL;
  char                 *ovlName     = NULL;
  char                 *outPrefix   = NULL;

  dumpParameters        params;

  gkRead_version        clearType = gkRead_latest;
  ovOverlapDisplayType  displType = ovOverlapAsCoords;
  char                  ovlString[1024];

  bool                  asOverlaps  = true;    //  What to show?
  bool                  asMetadata  = false;
  bool                  asCounts    = false;
  bool                  asErateLen  = false;

  bool                  asCoords = true;       //  How to show overlaps?
  bool                  asHangs  = false;
  bool                  asRaw    = false;
  bool                  asPAF    = false;
  bool                  asBinary = false;

#if 0
  if (asBinary) {
    char  binaryName[FILENAME_MAX];

    sprintf(binaryName, "%s.ovb", outPrefix);

    binaryFile = new ovFile(gkpStore, binaryName, ovFileFullWrite);
  }
#endif

  uint32                bgnID       = 1;
  uint32                endID       = UINT32_MAX;

  bool                  beVerbose   = false;


  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg=1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-G") == 0)
      gkpName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      ovlName = argv[++arg];


    else if (strcmp(argv[arg], "-overlaps") == 0) {
      asOverlaps = true;
      asMetadata = false;
      asCounts   = false;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if (strcmp(argv[arg], "-metadata") == 0) {
      asOverlaps = false;
      asMetadata = true;
      asCounts   = false;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if (strcmp(argv[arg], "-counts") == 0) {
      asOverlaps = false;
      asMetadata = false;
      asCounts   = true;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))        //  NOT CORRECT
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if (strcmp(argv[arg], "-eratelen") == 0) {
      asOverlaps = false;
      asMetadata = false;
      asCounts   = false;
      asErateLen = true;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))        //  NOT CORRECT
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }


    else if (strcmp(argv[arg], "-raw") == 0)
      clearType = gkRead_raw;

    else if (strcmp(argv[arg], "-obt") == 0)
      clearType = gkRead_corrected;

    else if (strcmp(argv[arg], "-utg") == 0)
      clearType = gkRead_trimmed;


    else if (strcmp(argv[arg], "-coords") == 0)
      asCoords = true;

    else if (strcmp(argv[arg], "-hangs") == 0)
      asHangs = true;

    else if (strcmp(argv[arg], "-raw") == 0)
      asRaw = true;

    else if (strcmp(argv[arg], "-paf") == 0)
      asPAF = true;

    else if (strcmp(argv[arg], "-binary") == 0)
      asBinary = true;


    else if (strcmp(argv[arg], "-no5p") == 0)
      params.no5p = true;

    else if (strcmp(argv[arg], "-no3p") == 0)
      params.no3p = true;

    else if (strcmp(argv[arg], "-nocontainer") == 0)
      params.noContainer = true;

    else if (strcmp(argv[arg], "-nocontained") == 0)
      params.noContained = true;

    else if (strcmp(argv[arg], "-noredundant") == 0)
      params.noRedundant = true;

    else if (strcmp(argv[arg], "-query") == 0)
      AS_UTL_decodeRange(argv[++arg], params.queryMin, params.queryMax);

    else if (strcmp(argv[arg], "-erate") == 0)
      AS_UTL_decodeRange(argv[++arg], params.erateMin, params.erateMax);

    else if (strcmp(argv[arg], "-length") == 0)
      AS_UTL_decodeRange(argv[++arg], params.lengthMin, params.lengthMax);


    else if (strcmp(argv[arg], "-v") == 0)
      beVerbose = true;

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (gkpName == NULL)
    err.push_back("ERROR: no input gkpStore (-G) supplied.\n");

  if (ovlName == NULL)
    err.push_back("ERROR: no input ovlStore (-O) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  gkStore   *gkpStore = gkStore::gkStore_open(gkpName);
  ovStore   *ovlStore = new ovStore(ovlName, gkpStore);
  uint32     ovlLen   = 0;
  uint32     ovlMax   = 65536;
  ovOverlap *ovl      = ovOverlap::allocateOverlaps(gkpStore, ovlMax);

  if (endID > gkpStore->gkStore_getNumReads())
    endID = gkpStore->gkStore_getNumReads();

  if (endID < bgnID)
    fprintf(stderr, "ERROR: invalid bgn/end range bgn=%u end=%u; only %u reads in the store\n", bgnID, endID, gkpStore->gkStore_getNumReads()), exit(1);

  ovlStore->setRange(bgnID, endID);



  //
  //  If dumping metadata, no filtering is needed, just tell the store to dump.
  //

  if (asMetadata) {
    ovlStore->dumpMetaData(bgnID, endID);
  }

  //
  //  If dumping as counts, maybe we need to filter, maybe not.
  //

  if ((asCounts) && (params.parametersAreDefaults() == true)) {
    uint32   *nopr = ovlStore->numOverlapsPerRead();

    for (uint32 ii=bgnID; ii<=endID; ii++)
      fprintf(stdout, "%u\t%u\n", ii, nopr[ii]);

    delete [] nopr;
  }

  if ((asCounts) && (params.parametersAreDefaults() == false)) {
    uint32   *nopr = new uint32 [gkpStore->gkStore_getNumReads() + 1];

    ovlLen = ovlStore->loadBlockOfOverlaps(ovl, ovlMax);

    while (ovlLen > 0) {
      for (uint32 oo=0; oo<ovlLen; oo++) {
        if (params.filterOverlap(ovl + oo) == true)
          continue;

        nopr[ovl[oo].a_iid]++;
        //nopr[ovl[oo].b_iid]++;
      }

      ovlLen = ovlStore->loadBlockOfOverlaps(ovl, ovlMax);
    }

    for (uint32 ii=bgnID; ii<=endID; ii++)
      fprintf(stdout, "%u\t%u\n", ii, nopr[ii]);

    delete [] nopr;
  }

  //
  //  If as erate-v-length, again, maybe, maybe not.
  //

  if ((asErateLen) && (params.parametersAreDefaults() == true)) {
    ovStoreHistogram  *hist = ovlStore->getHistogram();

    hist->dumpEvalueLength(stdout);

    delete hist;
  }

  if ((asErateLen) && (params.parametersAreDefaults() == true)) {
    ovStoreHistogram *hist = new ovStoreHistogram();

    ovlLen = ovlStore->loadBlockOfOverlaps(ovl, ovlMax);

    while (ovlLen > 0) {
      for (uint32 oo=0; oo<ovlLen; oo++) {
        if (params.filterOverlap(ovl + oo) == true)
          continue;

        hist->addOverlap(ovl + oo);
      }

      ovlLen = ovlStore->loadBlockOfOverlaps(ovl, ovlMax);
    }

    hist->dumpEvalueLength(stdout);

    delete hist;
  }

  //
  //  But if dumping actual overlaps, we've got to filter, and
  //  change the output format willy nilly.
  //
  
  if (asOverlaps) {
    ovlLen = ovlStore->loadBlockOfOverlaps(ovl, ovlMax);

    while (ovlLen > 0) {
      for (uint32 oo=0; oo<ovlLen; oo++) {
        if (params.filterOverlap(ovl + oo) == true)
          continue;

        if      (asCoords) {
          fputs(ovl[oo].toString(ovlString, ovOverlapAsCoords, true), stdout);
        }

        else if (asHangs) {
          fputs(ovl[oo].toString(ovlString, ovOverlapAsHangs, true), stdout);
        }

        else if (asRaw) {
          fputs(ovl[oo].toString(ovlString, ovOverlapAsRaw, true), stdout);
        }

        else if (asPAF) {
          fputs(ovl[oo].toString(ovlString, ovOverlapAsPaf, true), stdout);
        }

        else if (asBinary) {
          //binaryFile->writeOverlap(&overlap);
        }

        else {
        }
      }

      ovlLen = ovlStore->loadBlockOfOverlaps(ovl, ovlMax);
    }
  }

  //
  //  Phew, that was a lot of work.
  //

  delete ovlStore;
  gkpStore->gkStore_close();

  exit(0);
}
