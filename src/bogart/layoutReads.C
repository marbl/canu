
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_PlaceContains.H"

#include "AS_BAT_Outputs.H"

#include "AS_BAT_TigGraph.H"



ReadInfo         *RI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;



void
importTigsFromReadList(char const *prefix,
                       TigVector  &tigs,
                       char const *readListPath) {
  uint32        lineLen = 0;
  uint32        lineMax = 0;
  char         *line    = nullptr;
  splitToWords  words;

  FILE         *listFile = AS_UTL_openInputFile(readListPath);

  while (AS_UTL_readLine(line, lineLen, lineMax, listFile) == true) {
    for (uint32 ii=0; ii<lineLen; ii++)    //  Change commas to spaces
      if (line[ii] == ',')                 //  so we can split out words
        line[ii] = ' ';                    //  easily.

    words.split(line);

    //  Make a new tig.

    Unitig *tig = tigs.newUnitig();

    //  Add all the reads to it.

    for (uint32 rr=1; rr<words.numWords(); rr++) {
      uint32   fi   = 0;
      int32    fbgn = 0;
      int32    fend = 0;
      int32    flen = 0;
      bool     ffwd = true;

      //  Iterate over every letter in the name.  On the first digit, extract
      //  the read ID.  And set 'fwd' based on the last =+' or '-'.

      for (uint32 pp=0; words[rr][pp] != 0; pp++) {
        if ((fi == 0) && (isdigit(words[rr][pp]))) {
          fi   = strtouint32(&words[rr][pp]);
          flen = RI->readLength(fi);
        }

        if (words[rr][pp] == '-')
          ffwd = false;

        if (words[rr][pp] == '+')
          ffwd = true;
      }

      if (fi == 0)
        writeLog("failed to parse read ID from '%s'\n", words[rr]);
      if (flen == 0)
        writeLog("reference to non-existent read ID %u from '%s'\n", fi, words[rr]);
      assert(fi   != 0);
      assert(flen != 0);

      //  If this is the first read, we can immediately add it.

      if (rr == 1) {
        fbgn = 0;
        fend = flen;

        if (ffwd == true) {
          writeLog("place initial read %8u len %5u --> at %7u-%-7u\n", fi, flen, fbgn, fend);
          tig->addRead(ufNode(fi, fbgn, fend));
        } else {
          writeLog("place initial read %8u len %5u <-- at %7u-%-7u\n", fi, flen, fend, fbgn);
          tig->addRead(ufNode(fi, fend, fbgn));
        }

        continue;
      }

      //  Find the overlap to the last read in the tig.

      uint32       pi       = tig->ufpath.back().ident;
      bool         pfwd     = tig->ufpath.back().isForward();
      ufNode      &pnode    = tig->ufpath.back();

      uint32       povlLen  = 0;
      BAToverlap  *povl     = OC->getOverlaps(pi, povlLen);

      bool         wantFlip = (pfwd != ffwd);
      bool         placed   = false;

      for (uint32 poo=0; poo<povlLen; poo++) {    //  Iterator 'poo' for consistency with 'povl'.  Not because it's funny.  You'll have to trust me on that.
        if (povl[poo].b_iid != fi)           //  Overlap not to the 
          continue;                          //  read we want to place.

        if (povl[poo].flipped != wantFlip)   //  Overlap flippedness doesn't
          continue;                          //  agree with what we need.

        //  Must be the overlap we want.

        placed = true;

        //  Compute the position of the read.  Some nuance to point out here.
        //  If pnode is forward, then .bgn is the min coord, and a_hang _should_ be positive.
        //  If pnode is reverse, then .bgn is the max coord, and a_hang _should_ be negative.

        fbgn = pnode.position.bgn + povl[poo].a_hang;
        fend = pnode.position.bgn + povl[poo].a_hang + flen;

        if (pfwd ==  true)   assert(povl[poo].a_hang >= 0);
        if (pfwd == false)   assert(povl[poo].a_hang <= 0);

        //  However, testing with bogart layouts occasionally resulted in the
        //  next read being placed before the previous read, probably due to
        //  the 'optimization' baloney.

        if (fbgn < pnode.position.min()) {
          fbgn = pnode.position.min();
          fend = pnode.position.min() + flen;
        }

        //  Now just compute the position of the new read and plop it in the
        //  tig.

        if (ffwd == true) {
          writeLog("place         read %8u len %5u --> at %7u-%-7u\n", fi, flen, fbgn, fend);
          tig->addRead(ufNode(fi, fbgn, fend));
        } else {
          writeLog("place         read %8u len %5u <-- at %7u-%-7u\n", fi, flen, fend, fbgn);
          tig->addRead(ufNode(fi, fend, fbgn));
        }
      }

      if (placed == false) {
        writeLog("Failed to place '%s' next read '%s' (flipped %d) with overlaps:\n", words[rr-1], words[rr], wantFlip);

        for (uint32 poo=0; poo<povlLen; poo++)
          writeLog("  A %8u B %8u hangs %5d,%5d flip %d\n",
                   povl[poo].a_iid, povl[poo].b_iid, povl[poo].a_hang, povl[poo].b_hang, povl[poo].flipped);

        flushLog();
        assert(0);
      }
    }
  }

  delete [] line;

  AS_UTL_closeFile(listFile, readListPath);
}



int
main (int argc, char **argv) {
  char const  *seqStorePath   = NULL;
  char const  *ovlStorePath   = NULL;
  char const  *readListPath   = NULL;
  char const  *prefix         = NULL;

  double       erateGraph     = 1e-5; //0.075;
  double       deviationGraph = 6.0;
  double       erateMax       = 1e-5; //0.100;

  uint32       minReadLen     = 0;
  uint32       maxReadLen     = UINT32_MAX;

  uint64       ovlCacheMemory = UINT64_MAX;
  uint32       minOverlapLen  = 500;

  uint64       genomeSize     = 0;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-S") == 0) {
      seqStorePath = argv[++arg];
    }

    else if (strcmp(argv[arg], "-O") == 0) {
      ovlStorePath = argv[++arg];
    }

    else if (strcmp(argv[arg], "-R") == 0) {
      readListPath = argv[++arg];
    }

    else if (strcmp(argv[arg], "-o") == 0) {
      prefix = argv[++arg];
    }

    else if (strcmp(argv[arg], "-gs") == 0) {
      genomeSize = strtouint64(argv[++arg]);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (genomeSize   == 0)       err.push_back("Genome size (-gs option) must be supplied\n");
  if (seqStorePath == NULL)    err.push_back("No sequence store (-S option) supplied.\n");
  if (ovlStorePath == NULL)    err.push_back("No overlap store (-O option) supplied.\n");
  if (readListPath == NULL)    err.push_back("No list of reads (-R option) supplied.\n");
  if (prefix       == NULL)    err.push_back("No output prefix name (-o option) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqPath -O ovlPath -R readListPath -o outPrefix ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Mandatory Parameters:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqPath        Mandatory path to an existing seqStore.\n");
    fprintf(stderr, "  -O ovlPath        Mandatory path to an existing ovlStore.\n");
    fprintf(stderr, "  -R readListPath   Mandatory path to an existing ovlStore.\n");
    fprintf(stderr, "  -o outPrefix      Mandatory prefix for the output files.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //

  writeStatus("\n");
  writeStatus("==> LOADING READ AND OVERLAP INFORMATION.\n");
  writeStatus("\n");

  setLogFile(prefix, "loadInformation");

  RI = new ReadInfo(seqStorePath, prefix, minReadLen, maxReadLen);
  OC = new OverlapCache(ovlStorePath, prefix, max(erateMax, erateGraph), minOverlapLen, ovlCacheMemory, genomeSize);

  //

  writeStatus("\n");
  writeStatus("==> READING LAYOUTS AND BUILDING TIGS.\n");
  writeStatus("\n");

  setLogFile(prefix, "layoutTigs");

  TigVector  contigs(RI->numReads());  //  Both initial greedy tigs and final contigs

  importTigsFromReadList(prefix, contigs, readListPath);

  setLogFile(prefix, "layoutTigsOpt");
  contigs.optimizePositions(prefix, "layoutTigsOpt");
  reportTigs(contigs, prefix, "layoutTigsOpt", genomeSize);

  //

  writeStatus("\n");
  writeStatus("==> INSERTING CONTAINED READS.\n");
  writeStatus("\n");

  setLogFile(prefix, "placeContains");

  set<uint32>   placedReads;

  placeUnplacedUsingAllOverlaps(contigs,
                                deviationGraph,
                                erateGraph,
                                prefix,
                                placedReads);

  setLogFile(prefix, "placeContainsOpt");
  contigs.optimizePositions(prefix, "placeContainsOpt");
  reportTigs(contigs, prefix, "placeContainsOpt", genomeSize);

  //

  writeStatus("\n");
  writeStatus("==> GENERATE OUTPUTS.\n");
  writeStatus("\n");

  setLogFile(prefix, "generateOutputs");

  reportTigs(contigs, prefix, "final", genomeSize);
 
  //setParentAndHang(contigs);
  writeTigsToStore(contigs, prefix, "ctg", true);

  //

  delete OC;
  delete RI;

  writeStatus("\n");
  writeStatus("Bye.\n");

  return(0);
}

