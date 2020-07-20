
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
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_PlaceContains.H"

#include "AS_BAT_Outputs.H"

#include <unordered_map>
#include <unordered_set>

int32  MAX_SKIP_LIMIT    = 1;

ReadInfo         *RI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;

uint32 parseReadID(char *word, 
                   int32 &flen,
                   bool &ffwd) {

   //  Iterate over every letter in the name.  On the first digit, extract
   //  the read ID.  And set 'fwd' based on the last =+' or '-'.

   uint32 fi = 0;
   for (uint32 pp=0; word[pp] != 0; pp++) {
      if ((fi == 0) && (isdigit(word[pp]))) {
         fi   = strtouint32(&word[pp]);
         flen = RI->readLength(fi);
       }

       if (word[pp] == '-')
          ffwd = false;

       if (word[pp] == '+')
         ffwd = true;
   }
   if (fi == 0)
     writeLog("failed to parse read ID from '%s'\n", word);
   if (flen == 0)
     writeLog("reference to non-existent read ID %u from '%s'\n", fi, word);
   assert(fi   != 0);
   assert(flen != 0);
   
   return fi;
}

void processLine(char *line, uint32 lineLen) {
    for (uint32 ii=0; ii<lineLen; ii++)    //  Change commas to spaces
      if (line[ii] == ',')                 //  so we can split out words
        line[ii] = ' ';                    //  easily.
}

void
importTigsFromReadList(char const *prefix,
                       TigVector  &tigs,
                       char const *readListPath,
                       uint32 seed) {
  uint32        lineLen = 0;
  uint32        lineMax = 0;
  char         *line    = nullptr;
  splitToWords  words;

  srand(seed);
  std::unordered_map<uint32, uint32> readMap;
  std::unordered_set<uint32> usedReads;

  FILE         *listFile = AS_UTL_openInputFile(readListPath);

  // we make two passes through the reads
  // the first counts the reads to have information to randomize them
  // the second will actually build the tig
  //
  while (AS_UTL_readLine(line, lineLen, lineMax, listFile) == true) {
    processLine(line, lineLen);

    words.split(line);

    for (uint32 rr=1; rr<words.numWords(); rr++) {
      uint32   fi   = 0;
      int32    flen = 0;
      bool     ffwd = true;

      fi = parseReadID(words[rr], flen, ffwd);

      if (readMap.find(fi) == readMap.end()) {
         readMap[fi] = 0;
      }
      readMap[fi] += 1;
    }
  }
  AS_UTL_closeFile(listFile, readListPath);

  // OK go now build the tig
  listFile = AS_UTL_openInputFile(readListPath);
  while (AS_UTL_readLine(line, lineLen, lineMax, listFile) == true) {
    processLine(line, lineLen);

    words.split(line);

    //  Make a new tig.

    Unitig *tig = tigs.newUnitig();

    //  Add all the reads to it.
    bool isFirst = true;
    uint32 lastUsed = 0;
    uint32 skipped = 0;
    uint32 notSkipped = 0;

    for (uint32 rr=1; rr<words.numWords(); rr++) {
      uint32   fi   = 0;
      int32    fbgn = 0;
      int32    fend = 0;
      int32    flen = 0;
      bool     ffwd = true;

      fi = parseReadID(words[rr], flen, ffwd);

      // if we already used this read keep going
      if (usedReads.find(fi) != usedReads.end())
         continue;

      // generate a random number based on how many times we saw this read
      uint32 r = (rand() % readMap[fi]) + 1;
      //fprintf(stderr, "processing read %d with previous count %d and random number generated %d\n", fi, readMap[fi], r);
      // if we match randomly or we already skipped too many reads recently, we use this read
      if ((readMap[fi] > 1 && notSkipped > MAX_SKIP_LIMIT) || (r != 1 && skipped < MAX_SKIP_LIMIT)) {
         //fprintf(stderr, "Skpping read %d with %d count because value was %d\n", fi, readMap[fi], r);
         skipped++;
         notSkipped = 0;
         continue;
      }
      skipped =0 ;
      notSkipped++;
      usedReads.insert(fi);

      //  If this is the first read, we can immediately add it.
      if (isFirst) {
        lastUsed = rr;
        isFirst = false;
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

        if ((pfwd == true && povl[poo].a_hang < -1) || (pfwd == false && povl[poo].a_hang > 1))      //Overlap on the expected side of the read
           continue;

       //  Must be the overlap we want.

        placed = true;
        lastUsed = rr;

        //  Compute the position of the read.  Some nuance to point out here.
        //  If pnode is forward, then .bgn is the min coord, and a_hang _should_ be positive.
        //  If pnode is reverse, then .bgn is the max coord, and a_hang _should_ be negative.

        fbgn = pnode.position.bgn + povl[poo].a_hang;
        fend = pnode.position.bgn + povl[poo].a_hang + flen;

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
        writeLog("Failed to place '%s' next read '%s' (flipped %d) with overlaps:\n", words[lastUsed], words[rr], wantFlip);

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
  char const  *seqStorePath       = NULL;
  char const  *ovlStorePath       = NULL;
  char const  *readListPath       = NULL;
  char const  *prefix             = NULL;

  double       erateGraph         = 1e-5;
  double       deviationGraph     = 6.0;
  double       erateMax           = 1e-5;

  uint32       minReadLen         = 0;
  uint32       maxReadLen         = UINT32_MAX;

  uint64       ovlCacheMemory     = UINT64_MAX;
  uint32       minOverlapLen      = 500;

  uint64       genomeSize         = 0;
  bool         doContainPlacement = true;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  uint32 seed = time(NULL);

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

    else if (strcmp(argv[arg], "-eM") == 0) {
      erateMax = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-eg") == 0) {
      erateGraph = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-seed") == 0) {
      seed = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-nocontains") == 0) {
      doContainPlacement = false;
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
    fprintf(stderr, "  -S  seqPath        Mandatory path to an existing seqStore.\n");
    fprintf(stderr, "  -O  ovlPath        Mandatory path to an existing ovlStore.\n");
    fprintf(stderr, "  -R  readListPath   Mandatory path to an existing ovlStore.\n");
    fprintf(stderr, "  -gs genomeSize     Mandatory genome size in bp.\n");
    fprintf(stderr, "  -o  outPrefix      Mandatory prefix for the output files.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -eM erate          Max error rate of overlaps to load.\n");
    fprintf(stderr, "  -eg erate          Max error rate of overlaps to use for placing contained reads.\n");
    fprintf(stderr, "  -nocontains        Do not place contained reads.\n");

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

  importTigsFromReadList(prefix, contigs, readListPath, seed);

  setLogFile(prefix, "layoutTigsOpt");
  contigs.optimizePositions(prefix, "layoutTigsOpt");
  reportTigs(contigs, prefix, "layoutTigsOpt", genomeSize);

  //


  if (doContainPlacement == false) {
    writeStatus("\n");
    writeStatus("==> INSERTING CONTAINED READS DISABLED BY OPTION -nocontains.\n");
    writeStatus("\n");
  }

  else {
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
  }

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

