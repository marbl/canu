
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

#include "sqStore.H"
#include "ovStore.H"

#include <algorithm>
using namespace std;



class bogartStatus {
public:
  uint64  best5id      : 29;
  uint64  best53p      : 1;  // Unwieldy - best edge from my 5' is to the 3' of 'best5id'.

  uint64  best3id      : 29;
  uint64  best33p      : 1;

  uint64  unused       : 1;
  uint64  isSingleton  : 1;
  uint64  isContained  : 1;
  uint64  isSuspicious : 1;
};



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

    status             = NULL;


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

  ~dumpParameters() {
    delete [] status;
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

  void        loadBogartStatus(const char *prefix, uint32 nReads);

  void        drawPicture(uint32         Aid,
                          ovOverlap     *overlaps,
                          uint64         overlapsLen,
                          sqStore       *seqStore,
                          bool           withScores);


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

    if ((noBogartContained) && (status) && (status[overlap->b_iid].isContained == true)) {
      filtered = true;
    }

    if ((noBogartSuspicious) && (status) && (status[overlap->b_iid].isSuspicious == true)) {
      filtered = true;
    }

    if ((noBogartSingleton) && (status) && (status[overlap->b_iid].isSingleton == true)) {
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

  bool           no5p;
  bool           no3p;
  bool           noContainer;
  bool           noContained;
  bool           noRedundant;

  bool           noBogartContained;
  bool           noBogartSuspicious;
  bool           noBogartSingleton;

  uint32         queryMin;
  uint32         queryMax;

  double         erateMin;
  double         erateMax;

  uint32         lengthMin;
  uint32         lengthMax;

  char          *bogartPath;
  bogartStatus  *status;

  //  Counts of what we filtered.

  uint64         ovlKept;
  uint64         ovlFiltered;

  uint64         ovl5p;
  uint64         ovl3p;
  uint64         ovlContainer;
  uint64         ovlContained;
  uint64         ovlRedundant;

  uint64         ovlErateHi;
  uint64         ovlErateLo;

  uint64         ovlLengthHi;
  uint64         ovlLengthLo;
};



void
dumpParameters::loadBogartStatus(const char *prefix, uint32 nReads) {
  char          EN[FILENAME_MAX+1];
  char          SN[FILENAME_MAX+1];
  char          GN[FILENAME_MAX+1];
  char          NN[FILENAME_MAX+1];
  splitToWords  W;

  if (prefix == NULL)
    return;

  snprintf(EN, FILENAME_MAX, "%s.edges", prefix);
  snprintf(SN, FILENAME_MAX, "%s.edges.suspicious", prefix);
  snprintf(GN, FILENAME_MAX, "%s.singletons", prefix);

  allocateArray(status, nReads+1, resizeArray_clearNew);

  FILE *E = AS_UTL_openInputFile(EN);
  fgets(NN, FILENAME_MAX, E);
  while (!feof(E)) {
    W.split(NN);

    uint32  id = W(0);

    status[id].best5id = W(2);
    status[id].best53p = (W[3][0] == '3');

    status[id].best3id = W(4);
    status[id].best33p = (W[5][0] == '3');

    status[id].isSingleton  = false;
    status[id].isContained  = ((W.numWords() > 10) && (W[10][0] == 'c'));
    status[id].isSuspicious = false;

    fgets(NN, FILENAME_MAX, E);
  }
  AS_UTL_closeFile(E, EN);


  FILE *S = AS_UTL_openInputFile(SN);
  fgets(NN, FILENAME_MAX, S);
  while (!feof(S)) {
    W.split(NN);

    uint32  id = W(0);

    status[id].best5id = W(2);
    status[id].best53p = (W[3][0] == '3');

    status[id].best3id = W(4);
    status[id].best33p = (W[5][0] == '3');

    status[id].isSingleton  = false;
    status[id].isContained  = ((W.numWords() > 10) && (W[10][0] == 'c'));
    status[id].isSuspicious = true;

    fgets(NN, FILENAME_MAX, S);
  }
  AS_UTL_closeFile(S, SN);


  FILE *G = AS_UTL_openInputFile(GN);
  fgets(NN, FILENAME_MAX, G);
  while (!feof(G)) {
    W.split(NN);

    uint32  id = W(0);

    status[id].best5id = 0;
    status[id].best53p = 0;

    status[id].best3id = 0;
    status[id].best33p = 0;

    status[id].isSingleton  = true;
    status[id].isContained  = false;
    status[id].isSuspicious = false;

    fgets(NN, FILENAME_MAX, G);
  }
  AS_UTL_closeFile(G, GN);
}


class sortByPosition {
public:
  bool operator()(const ovOverlap &a, const ovOverlap &b) {
    if (a.a_bgn() < b.a_bgn())  return(true);
    if (a.a_bgn() > b.a_bgn())  return(false);

    return(a.a_end() < b.a_end());
  };
};



void
dumpParameters::drawPicture(uint32         Aid,
                            ovOverlap     *overlaps,
                            uint64         overlapsLen,
                            sqStore       *seqStore,
                            bool           withScores) {
  char     line[256] = {0};

  uint32   MHS   = 9;  //  Max Hang Size, amount of padding for "+### "

  sqRead   *A    = seqStore->sqStore_getRead(Aid);
  uint32    Alen = A->sqRead_sequenceLength();

  for (int32 i=0; i<256; i++)
    line[i] = ' ';

  for (int32 i=0; i<100; i++)
    line[i + MHS] = '-';
  line[ 99 + MHS] = '>';
  line[100 + MHS] = 0;

  fprintf(stdout, "A %7d:%-7d A %9d %7d:%-7d %7d          %s %s%s\n",
          0, Alen,
          Aid,
          0, Alen, Alen,
          line,
          ((status) && (status[Aid].isContained))  ? "contained"  : "",
          ((status) && (status[Aid].isSuspicious)) ? "suspicious" : "");

  sort(overlaps, overlaps + overlapsLen, sortByPosition());

  //  Build ascii representations for each overlapping read.

  for (uint32 o=0; o<overlapsLen; o++) {
    uint32    Bid  = overlaps[o].b_iid;
    sqRead   *B    = seqStore->sqStore_getRead(Bid);
    uint32    Blen = B->sqRead_sequenceLength();

    //  Find bgn/end points on each read.  If the overlap is reverse complement,
    //  the B coords are flipped so that bgn > end.

    uint32   ovlBgnA = overlaps[o].a_bgn();
    uint32   ovlEndA = overlaps[o].a_end();

    uint32   ovlBgnB = overlaps[o].b_bgn();
    uint32   ovlEndB = overlaps[o].b_end();

    assert(ovlBgnA < ovlEndA);  //  The A coordiantes are always forward

    if (overlaps[o].flipped() == false)
      assert(ovlBgnB < ovlEndB);  //  Forward overlaps are forward
    else
      assert(ovlEndB < ovlBgnB);  //  Flipped overlaps are reversed

    //  For the A read, find the points in our string representation where the overlap ends.

    uint32 ovlStrBgn = (int32)floor(ovlBgnA * 100.0 / Alen + MHS);
    uint32 ovlStrEnd = (int32)ceil (ovlEndA * 100.0 / Alen + MHS);

    //  Fill the string representation with spaces, then fill the string with dashes where the read
    //  is, add an arrow, and terminate the string.

    for (int32 i=0; i<256; i++)
      line[i] = ' ';

    //  Decide how to draw this overlap.
    //    For best edges, use '='.
    //    For contained,  use '-', alternating with spaces.
    //    For suspicious, use '*', alternating with dashes.
    //    For edges,      use '-', solid.

    bool  isBest = ((status) && (((status[Aid].best5id == Bid) && (overlaps[o].overlapAEndIs5prime() == true)) ||
                                 ((status[Aid].best3id == Bid) && (overlaps[o].overlapAEndIs3prime() == true))));
    bool  isCont = ((status) && (status[Bid].isContained));
    bool  isSusp = ((status) && (status[Bid].isSuspicious));

    //  This bit of confusion makes sure that the alternating overlap lines (especially '- - - -')
    //  end with a dash.

    bool  oddEven = (overlaps[o].flipped() == false) ? (false) : (((ovlStrEnd - ovlStrBgn) % 2) == false);

    if      (isCont == true) {
      for (uint32 i=ovlStrBgn; i<ovlStrEnd; i++)
        line[i] = (oddEven = !oddEven) ? '-' : ' ';
    }

    else if (isSusp == true) {
      for (uint32 i=ovlStrBgn; i<ovlStrEnd; i++)
        line[i] = (oddEven = !oddEven) ? '-' : '*';
    }

    else {
      char  c = (isBest) ? '=' : '-';

      for (uint32 i=ovlStrBgn; i<ovlStrEnd; i++)
        line[i] = c;
    }

    if (overlaps[o].flipped() == true)
      line[ovlStrBgn] = '<';
    else
      line[ovlStrEnd-1] = '>';

    assert(line[ovlStrBgn]   != ' ');
    assert(line[ovlStrEnd-1] != ' ');

    line[ovlStrEnd] = 0;

    //  For the B read, find how much is unaliged on each end.  Though the store directly keeps this information,
    //  we can't get to it, and have to reverse the compuitation.

    uint32  ovlBgnHang = 0;
    uint32  ovlEndHang = 0;

    if (overlaps[o].flipped() == false) {
      ovlBgnHang = ovlBgnB;
      ovlEndHang = Blen - ovlEndB;
    } else {
      ovlBgnHang = Blen - ovlBgnB;
      ovlEndHang = ovlEndB;
    }

    //  Paste the bgn hang into the overlap string.
    if (ovlBgnHang > 0) {
      char  str[256];
      int32 len;

      snprintf(str, 256, "+%d", ovlBgnHang);
      len = strlen(str);

      for (int32 i=0; i<len; i++)
        line[ovlStrBgn - len - 1 + i] = str[i];
    }

    //  Append the end hang.
    if (ovlEndHang > 0) {
      snprintf(line + ovlStrEnd, 256 - ovlStrEnd, " +%d", ovlEndHang);
    }

    //  Set flags for best edge and singleton/contained/suspicious.

    char  olapClass[4] = { 0 };

    if ((status) && (status[Aid].best5id == Bid) && (overlaps[o].overlapAEndIs5prime() == true)) {
      olapClass[0] = ' ';
      olapClass[1] = ' ';
      olapClass[2] = 'B';
    }

    if ((status) && (status[Aid].best3id == Bid) && (overlaps[o].overlapAEndIs3prime() == true)) {
      olapClass[0] = ' ';
      olapClass[1] = ' ';
      olapClass[2] = 'B';
    }

    if (olapClass[2] == 'B')
      for (uint32 ii=0; line[ii]; ii++)
        if (line[ii] == '-')
          line[ii] = '=';

    if ((status) && (status[Bid].isSingleton)) {
      olapClass[0] = ' ';
      olapClass[1] = 'S';
    }

    if ((status) && (status[Bid].isContained)) {
      olapClass[0] = ' ';
      olapClass[1] = 'C';
    }

    if ((status) && (status[Bid].isSuspicious)) {
      olapClass[0] = ' ';
      olapClass[1] = '!';
    }

    //  Report!

    if (withScores)
      fprintf(stdout, "A %7d:%-7d B %9d %7d:%-7d %7d %7hu %5.2f%%  %s%s\n",
              ovlBgnA,
              ovlEndA,
              Bid,
              min(ovlBgnB, ovlEndB),
              max(ovlBgnB, ovlEndB),
              Blen,
              overlaps[o].overlapScore(),
              overlaps[o].erate() * 100.0,
              line,
              olapClass);
    else
      fprintf(stdout, "A %7d:%-7d B %9d %7d:%-7d %7d  %5.2f%%  %s%s\n",
              ovlBgnA,
              ovlEndA,
              Bid,
              min(ovlBgnB, ovlEndB),
              max(ovlBgnB, ovlEndB),
              Blen,
              overlaps[o].erate() * 100.0,
              line,
              olapClass);
  }
}



int
main(int argc, char **argv) {
  char                 *seqName     = NULL;
  char                 *ovlName     = NULL;
  char                 *outPrefix   = NULL;
  char                 *bogartPath  = NULL;

  dumpParameters        params;

  char                  ovlString[1024];

  bool                  asOverlaps  = true;    //  What to show?
  bool                  asPicture   = false;
  bool                  asMetadata  = false;
  bool                  asCounts    = false;
  bool                  asErateLen  = false;

  bool                  asCoords    = true;       //  How to show overlaps?
  bool                  asHangs     = false;
  bool                  asRaw       = false;
  bool                  asPAF       = false;
  bool                  asBinary    = false;
  bool                  withScores  = false;

#if 0
  if (asBinary) {
    char  binaryName[FILENAME_MAX];

    sprintf(binaryName, "%s.ovb", outPrefix);

    binaryFile = new ovFile(seqStore, binaryName, ovFileFullWrite);
  }
#endif

  uint32                bgnID       = 1;
  uint32                endID       = UINT32_MAX;

  bool                  beVerbose   = false;


  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg=1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-S") == 0)
      seqName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      ovlName = argv[++arg];


    else if (strcmp(argv[arg], "-overlaps") == 0) {
      asOverlaps = true;
      asPicture  = false;
      asMetadata = false;
      asCounts   = false;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if ((strcmp(argv[arg], "-picture") == 0) ||
             (strcmp(argv[arg], "-p")       == 0)) {
      asOverlaps = false;
      asPicture  = true;
      asMetadata = false;
      asCounts   = false;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if (strcmp(argv[arg], "-metadata") == 0) {
      asOverlaps = false;
      asPicture  = false;
      asMetadata = true;
      asCounts   = false;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if (strcmp(argv[arg], "-counts") == 0) {
      asOverlaps = false;
      asPicture  = false;
      asMetadata = false;
      asCounts   = true;
      asErateLen = false;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))        //  NOT CORRECT
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }

    else if (strcmp(argv[arg], "-eratelen") == 0) {
      asOverlaps = false;
      asPicture  = false;
      asMetadata = false;
      asCounts   = false;
      asErateLen = true;
      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))        //  NOT CORRECT
        AS_UTL_decodeRange(argv[++arg], bgnID, endID);
    }


    else if (strcmp(argv[arg], "-raw") == 0)
      sqRead_setDefaultVersion(sqRead_raw);

    else if (strcmp(argv[arg], "-obt") == 0)
      sqRead_setDefaultVersion(sqRead_corrected);

    else if (strcmp(argv[arg], "-utg") == 0)
      sqRead_setDefaultVersion(sqRead_trimmed);


    else if (strcmp(argv[arg], "-coords") == 0) {
      asCoords = true;
      asHangs  = false;
      asRaw    = false;
      asPAF    = false;
      asBinary = false;
    }

    else if (strcmp(argv[arg], "-hangs") == 0) {
      asCoords = false;
      asHangs  = true;
      asRaw    = false;
      asPAF    = false;
      asBinary = false;
    }

    else if (strcmp(argv[arg], "-raw") == 0) {
      asCoords = false;
      asHangs  = false;
      asRaw    = true;
      asPAF    = false;
      asBinary = false;
    }

    else if (strcmp(argv[arg], "-paf") == 0) {
      asCoords = false;
      asHangs  = false;
      asRaw    = false;
      asPAF    = true;
      asBinary = false;
    }

    else if (strcmp(argv[arg], "-binary") == 0) {
      asCoords = false;
      asHangs  = false;
      asRaw    = false;
      asPAF    = false;
      asBinary = true;
    }


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

    else if (strcmp(argv[arg], "-bogart") == 0)
      bogartPath = argv[++arg];


    else if (strcmp(argv[arg], "-v") == 0)
      beVerbose = true;

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqName == NULL)
    err.push_back("ERROR: no input seqStore (-S) supplied.\n");

  if (ovlName == NULL)
    err.push_back("ERROR: no input ovlStore (-O) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //
  //  Open stores, allocate space to store overlaps.
  //

  sqStore   *seqStore = sqStore::sqStore_open(seqName);
  ovStore   *ovlStore = new ovStore(ovlName, seqStore);

  uint32     ovlLen   = 0;
  uint32     ovlMax   = 65536;
  ovOverlap *ovl      = ovOverlap::allocateOverlaps(seqStore, ovlMax);

  //
  //  Fix up ranges and restrict the overlaps.
  //

  if (endID > seqStore->sqStore_getNumReads())
    endID = seqStore->sqStore_getNumReads();

  if (endID < bgnID)
    fprintf(stderr, "ERROR: invalid bgn/end range bgn=%u end=%u; only %u reads in the store\n", bgnID, endID, seqStore->sqStore_getNumReads()), exit(1);

  ovlStore->setRange(bgnID, endID);

  //
  //  Load bogart status.
  //

  params.loadBogartStatus(bogartPath, seqStore->sqStore_getNumReads());


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
    uint32   *nopr = new uint32 [seqStore->sqStore_getNumReads() + 1];

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
    ovStoreHistogram *hist = new ovStoreHistogram(seqStore);

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
  //  If as a picture, we need to build a list of the overlaps to show,
  //  then draw the picture.
  //

  if (asPicture) {
    for (uint32 rr=bgnID; rr<=endID; rr++) {
      uint32  ovlSav = 0;

      ovlLen = ovlStore->loadOverlapsForRead(rr, ovl, ovlMax);

      for (uint32 oo=0; oo<ovlLen; oo++)
        if (params.filterOverlap(ovl + oo) == false)   //  If not filtered,
          ovl[ovlSav++] = ovl[oo];                     //  save the overlap for drawing

      if (ovlSav > 0)
        params.drawPicture(rr, ovl, ovlSav, seqStore, withScores);
    }
  }

  //
  //  Phew, that was a lot of work.
  //

  delete ovlStore;
  seqStore->sqStore_close();

  exit(0);
}
