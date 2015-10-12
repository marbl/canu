
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
 *    src/AS_CNS/tigStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-OCT-05 to 2014-MAR-31
 *      are Copyright 2009-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Michael Schatz from 2010-FEB-27 to 2010-MAY-17
 *      are Copyright 2010 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren beginning on 2014-APR-13
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2014-OCT-09
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

const char *mainid = "$Id$";

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

#include "AS_UTL_decodeRange.H"
#include "intervalList.H"

#include "tgTigSizeAnalysis.H"



#define DUMP_PROPERTIES       1
#define DUMP_CONSENSUS        4
#define DUMP_CONSENSUSGAPPED  5
#define DUMP_LAYOUT           6
#define DUMP_MULTIALIGN       7
#define DUMP_SIZES            8
#define DUMP_COVERAGE         9
#define DUMP_THINOVERLAP     10
#define DUMP_FMAP            11

#define OPERATION_LIST        1
#define OPERATION_PROPERTIES  2
#define OPERATION_TIG         3
#define OPERATION_EDIT        4
#define OPERATION_REPLACE     5
#define OPERATION_BUILD       6
#define OPERATION_COMPRESS    7




void
changeProperties(tgStore *tigStore,
                 char    *editName) {

#warning changeProperties not implemented
  fprintf(stderr, "changeProperties() not implemented\n");
  exit(1);

#if 0
  char  editLine[1024];

  errno = 0;
  FILE *editFile = fopen(editName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", editName, strerror(errno)), exit(1);

  fgets(editLine, 1024, editFile);
  while (!feof(editFile)) {
    chomp(editLine);

    //  Decode the string into three pieces, the operation, the tig id and the value.  Note that
    //  these are just pointers to the start of the piece, and that the pieces are NOT
    //  zero-terminated (so we can print an error of the whole line if needed).
    //
    char   *op = editLine;
    char   *tp = editLine;
    char   *vp = editLine;

    //  Skip whitespace before the operation
    while (isspace(*op))
      op++;

    //  Skip the operation, then skip whitespace before the tig.
    tp = op;
    while (!isspace(*tp))
      tp++;
    while (isspace(*tp))
      tp++;

    //  Skip the tig, then skip whitespace before the value.
    vp = tp;
    while (!isspace(*vp))
      vp++;
    while (isspace(*vp))
      vp++;

    int32   tid = atoi(tp);

    if        (strncmp(op, "unitig_coverage_stat", 20) == 0) {
      tigStore->setUnitigCoverageStat(tid, atof(vp));

    //} else if (strncmp(op, "unitig_microhet_prob", 20) == 0) {
    //tigStore->setUnitigMicroHetProb(tid, atof(vp));

    } else if (strncmp(op, "unitig_status", 13) == 0) {
      UnitigStatus  st = tigStore->getUnitigStatus(tid);
      switch (*vp) {
        case AS_UNIQUE:      st = AS_UNIQUE;      break;
        case AS_NOTREZ:      st = AS_NOTREZ;      break;
        case AS_SEP:         st = AS_SEP;         break;
        case AS_UNASSIGNED:  st = AS_UNASSIGNED;  break;
        default:
          fprintf(stderr, "unknown unitig_status in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigStatus(tid, st);

    } else if (strncmp(op, "unitig_suggest_repeat", 21) == 0) {
      bool ur = tigStore->getUnitigSuggestRepeat(tid);
      switch (*vp) {
        case 'T':  ur = true;    break;
        case 'F':  ur = false;   break;
        default:
          fprintf(stderr, "unknown unitig_suggest_repeat in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigSuggestRepeat(tid, ur);

    } else if (strncmp(op, "unitig_suggest_unique", 21) == 0) {
      bool ur = tigStore->getUnitigSuggestUnique(tid);
      switch (*vp) {
        case 'T':  ur = true;    break;
        case 'F':  ur = false;   break;
        default:
          fprintf(stderr, "unknown unitig_suggest_unique in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigSuggestUnique(tid, ur);

    } else if (strncmp(op, "unitig_force_repeat", 21) == 0) {
      bool ur = tigStore->getUnitigForceRepeat(tid);
      switch (*vp) {
        case 'T':  ur = true;    break;
        case 'F':  ur = false;   break;
        default:
          fprintf(stderr, "unknown unitig_force_repeat in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigForceRepeat(tid, ur);

    } else if (strncmp(op, "unitig_force_unique", 21) == 0) {
      bool ur = tigStore->getUnitigForceUnique(tid);
      switch (*vp) {
        case 'T':  ur = true;    break;
        case 'F':  ur = false;   break;
        default:
          fprintf(stderr, "unknown unitig_force_unique in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigForceUnique(tid, ur);

    } else if (strncmp(op, "contig_status", 13) == 0) {
      ContigStatus  st = tigStore->getContigStatus(tid);
      switch (*vp) {
        case AS_PLACED:      st = AS_PLACED;    break;
        case AS_UNPLACED:    st = AS_UNPLACED;  break;
        default:
          fprintf(stderr, "unknown contig_status in '%s'\n", editLine);
          break;
      }
      tigStore->setContigStatus(tid, st);

    } else {
      fprintf(stderr, "unknown edit '%s'\n", editLine);
    }

    fgets(editLine, 1024, editFile);
  }
#endif
}



void
dumpProperties(tgStore *tigStore,
               tgTig   *tig) {

  fprintf(stdout, "tigID            "F_U32"\n", tig->_tigID);
  fprintf(stdout, "coverageStat     %f\n",      tig->_coverageStat);
  //fprintf(stdout, "microhetProb     %f\n",      tig->_microhetProb);
  fprintf(stdout, "suggestRepeat    %d\n",      tig->_suggestRepeat);
  fprintf(stdout, "suggestUnique    %d\n",      tig->_suggestUnique);
  fprintf(stdout, "suggestCircular  %d\n",      tig->_suggestCircular);
  fprintf(stdout, "suggestHaploid   %d\n",      tig->_suggestHaploid);
  fprintf(stdout, "numChildren      "F_U32"\n", tig->_childrenLen);

#if GCCONTENT
  float gcContent = 0.0;
  int ulen = 0;
  int glen = 0;

  if (tig->consensus) {
    char *cns = Getchar(tig->consensus, 0);

    if (cns && *cns) {
      int gcs = 0;

      while (*cns) {
        glen++;

        if (*cns != '-') {
          if (*cns == 'G' || *cns == 'C') { gcs++; }
          ulen++;
        }

        cns++;
      }

      if (ulen) { gcContent = (float) gcs / (float) ulen; }
    }
  }

  fprintf(stdout, "gcContent           %f\n",      gcContent);
  fprintf(stdout, "uLen                %d\n",      ulen);
  fprintf(stdout, "gLen                %d\n",      glen);
#endif

  //  This dumped store-private info, like present, deleted, partition, version and file offset
  //tigStore->dumpMultiAlignR(tig->tigID());
}



void
dumpConsensus(tgStore *tigStore,
              tgTig   *tig,
              bool     withGaps,
              uint32   minCoverage) {

  if (tig->gappedLength() == 0)
    return;

  char   *cns = tig->gappedBases();
  uint32  len = tig->gappedLength();

  //  If a minCoverage is specified, convert the low coverage bases to underscores, which will be
  //  filtered later.

  if (minCoverage > 0) {
    intervalList<int32>  allL;

    for (uint32 i=0; i<tig->numberOfChildren(); i++) {
      tgPosition *imp = tig->getChild(i);

      uint32      bgn = imp->min();
      uint32      end = imp->max();

      allL.add(bgn, end - bgn);
    }

    intervalList<int32>  ID(allL);

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
      if (ID.depth(ii) >= minCoverage)
        continue;

      for (uint32 pp=ID.lo(ii); pp<ID.hi(ii); pp++)
        cns[pp] = '_';
    }
  }

  //  Now filter out gaps in the consensus.

  if (withGaps == false) {
    char *o = cns;
    char *n = cns;

    len = 0;

    while (*n) {
      if (*n != '-') {
        *o++ = *n;
        len++;
      }
      n++;
    }

    *o = 0;
  }

  //  If no min coverage, we can just dump the consensus and be done.  Plus we output a few bits of
  //  useful info with the sequence.

  if (minCoverage == 0) {
    fprintf(stdout, ">tig"F_U32" len="F_U32" reads="F_U32" covStat=%.2f\n%s\n",
            tig->tigID(),
            len,
            tig->numberOfChildren(),
            //tig->microhetProb(),
            tig->coverageStat(),
            cns);
    return;
  }

  //  Otherwise, we need to find subsequences in the consensus.  The useful bits of info aren't
  //  valid anymore.

  uint32  subsequence   = 0;

  for (uint32 bgn=0; bgn<len; bgn++)
    if (cns[bgn] == '_')
      cns[bgn] = 0;

  for (uint32 bgn=0; bgn<len; bgn++) {
    while ((cns[bgn] == 0) && (bgn < len))
      bgn++;

    if (bgn >= len)
      break;

    uint32 end = bgn + 1;

    while ((cns[end] != 0) && (end < len))
      end++;

    fprintf(stdout, ">tig%d.%u bgn=%u end=%u len=%u\n%s\n",
            tig->tigID(), subsequence, bgn, end, end-bgn,
            cns + bgn);

    bgn = end + 1;

    subsequence++;
  }
}



void
dumpCoverage(tgStore  *tigStore,
             tgTig    *tig,
             uint32    minCoverage,
             uint32    maxCoverage,
             uint64   *coverageHistogram,
             uint32    coverageHistogramLen,
             char     *outPrefix) {
  intervalList<int32>  allL;

  uint32        maxPos = tig->layoutLength();

  for (uint32 i=0; i<tig->numberOfChildren(); i++) {
    tgPosition *imp = tig->getChild(i);

    int32   bgn = imp->min();
    int32   end = imp->max();

    allL.add(bgn, end - bgn);
  }

  intervalList<int32>   ID(allL);

  intervalList<int32>   minL;
  intervalList<int32>   maxL;

  uint32  maxDepth    = 0;
  double  aveDepth    = 0;
  double  sdeDepth    = 0;

#if 0
  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
    if ((ID.depth(ii) < minCoverage) && (ID.lo(ii) != 0) && (ID.hi(ii) != maxPos)) {
      fprintf(stderr, "tig %d low coverage interval %ld %ld max %u coverage %u\n",
              tig->tigID(), ID.lo(ii), ID.hi(ii), maxPos, ID.depth(ii));
      minL.add(ID.lo(ii), ID.hi(ii) - ID.lo(ii) + 1);
    }

    if (maxCoverage <= ID.depth(ii)) {
      fprintf(stderr, "tig %d high coverage interval %ld %ld max %u coverage %u\n",
              tig->tigID(), ID.lo(ii), ID.hi(ii), maxPos, ID.depth(ii));
      maxL.add(ID.lo(ii), ID.hi(ii) - ID.lo(ii) + 1);
    }
  }
#endif

  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
    if (ID.depth(ii) > maxDepth)
      maxDepth = ID.depth(ii);

    aveDepth += (ID.hi(ii) - ID.lo(ii) + 1) * ID.depth(ii);

    if (ID.depth(ii) < coverageHistogramLen)
      coverageHistogram[ID.depth(ii)] += ID.hi(ii) - ID.lo(ii) + 1;
    else
      fprintf(stderr, "deep coverage %d\n", ID.depth(ii));
  }

  aveDepth /= maxPos;

  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
    sdeDepth += (ID.hi(ii) - ID.lo(ii) + 1) * (ID.depth(ii) - aveDepth) * (ID.depth(ii) - aveDepth);
  }

  sdeDepth = sqrt(sdeDepth / maxPos);

  if (maxDepth > 1000)
    fprintf(stderr, "DEEP unitig %u of length %u with maxDepth %u\n",
            tig->tigID(), maxPos, maxDepth);

  allL.merge();
  minL.merge();
  maxL.merge();


#if 0
  if      ((minL.numberOfIntervals() > 0) && (maxL.numberOfIntervals() > 0))
    fprintf(stderr, "tig %d has %u intervals, %u regions below %u coverage and %u regions at or above %u coverage\n",
            tig->tigID(),
            allL.numberOfIntervals(),
            minL.numberOfIntervals(), minCoverage,
            maxL.numberOfIntervals(), maxCoverage);
  else if (minL.numberOfIntervals() > 0)
    fprintf(stderr, "tig %d has %u intervals, %u regions below %u coverage\n",
            tig->tigID(),
            allL.numberOfIntervals(),
            minL.numberOfIntervals(), minCoverage);
  else if (maxL.numberOfIntervals() > 0)
    fprintf(stderr, "tig %d has %u intervals, %u regions at or above %u coverage\n",
            tig->tigID(),
            allL.numberOfIntervals(),
            maxL.numberOfIntervals(), maxCoverage);
  else
    fprintf(stderr, "tig %d has %u intervals\n",
            tig->tigID(),
            allL.numberOfIntervals());
#endif

  if (outPrefix) {
    char  outName[FILENAME_MAX];

    sprintf(outName, "%s.tig%08u.depth", outPrefix, tig->tigID());

    FILE *outFile = fopen(outName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
      fprintf(outFile, "%d\t%u\n", ID.lo(ii),     ID.depth(ii));
      fprintf(outFile, "%d\t%u\n", ID.hi(ii) - 1, ID.depth(ii));
    }

    fclose(outFile);

    FILE *gnuPlot = popen("gnuplot > /dev/null 2>&1", "w");

    if (gnuPlot) {
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "set output '%s.tig%08u.png'\n", outPrefix, tig->tigID());
      fprintf(gnuPlot, "set xlabel 'position'\n");
      fprintf(gnuPlot, "set ylabel 'coverage'\n");
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "plot '%s.tig%08u.depth' using 1:2 with lines title 'tig %u length %u', \\\n",
              outPrefix,
              tig->tigID(),
              tig->tigID(), maxPos);
      fprintf(gnuPlot, "     %f title 'mean %.2f +- %.2f', \\\n", aveDepth, aveDepth, sdeDepth);
      fprintf(gnuPlot, "     %f title '' lt 0 lc 2, \\\n", aveDepth - sdeDepth);
      fprintf(gnuPlot, "     %f title '' lt 0 lc 2\n",     aveDepth + sdeDepth);

      fclose(gnuPlot);
    }
  }
}


void
dumpThinOverlap(tgStore *tigStore,
                tgTig   *tig,
                uint32   minOverlap) {
  intervalList<int32>  allL;
  intervalList<int32>  ovlL;

  uint32        maxPos = tig->layoutLength();

  for (uint32 i=0; i<tig->numberOfChildren(); i++) {
    tgPosition *imp = tig->getChild(i);

    int32   bgn = imp->min();
    int32   end = imp->max();

    allL.add(bgn, end - bgn);
    ovlL.add(bgn, end - bgn);
  }

  allL.merge();
  ovlL.merge(minOverlap);

  if (ovlL.numberOfIntervals() == 1)
    return;


  if (maxPos < 10000)
    return;

#if 1
  for (uint32 i=0; i<ovlL.numberOfIntervals(); i++)
    fprintf(stderr, "tig %u IL %d %d\n",
            tig->tigID(),
            ovlL.lo(i), ovlL.hi(i));
#endif

  intervalList<int32> badL;

  for (uint32 i=1; i<ovlL.numberOfIntervals(); i++) {
    assert(ovlL.lo(i) < ovlL.hi(i-1));

    badL.add(ovlL.lo(i), ovlL.hi(i-1) - ovlL.lo(i));
    //badL.add(ovlL.hi(i-1), ovlL.lo(i) - ovlL.hi(i-1));
  }

#if 1
  for (uint32 i=0; i<badL.numberOfIntervals(); i++)
    fprintf(stderr, "tig %u BAD %d %d\n",
            tig->tigID(),
            badL.lo(i), badL.hi(i));
#endif

  for (uint32 i=0; i<tig->numberOfChildren(); i++) {
    tgPosition *imp = tig->getChild(i);

    uint32   bgn = imp->min();
    uint32   end = imp->max();

    bool    report = false;

    for (uint32 oo=0; oo<badL.numberOfIntervals(); oo++)
      if ((badL.lo(oo) <= end) &&
          (bgn         <= badL.hi(oo))) {
        report = true;
        break;
      }

    if (report)
      fprintf(stderr, "tig %d frag %u %u-%u\n",
              tig->tigID(),
              imp->ident(), imp->bgn(), imp->end());
  }

  fprintf(stderr, "tig %d max %u has %u intervals, %u enforcing minimum overlap of %u\n",
          tig->tigID(), maxPos,
          allL.numberOfIntervals(),
          ovlL.numberOfIntervals(), minOverlap);
}





void
operationCompress(char *tigName, int tigVers) {
  tgStore    *tigStore  = new tgStore(tigName, tigVers);
  uint32      nErrors   = 0;
  uint32      nCompress = 0;

  //  Pass 0:  Fail if this isn't the latest version.  If we try to compress something that isn't the
  //  latest version, versions after this still point to the uncompressed tigs.

  //if (tigStore->

  //  Pass 1:  Check that we aren't going to pull a tig out of the future and place it in the past.

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (tigStore->getVersion(ti) > tigVers) {
      fprintf(stderr, "WARNING:  Attempt to move future unitig "F_U32" from version "F_U32" to previous version %d.\n",
              ti, tigStore->getVersion(ti), tigVers);
      nErrors++;
    } else if (tigStore->getVersion(ti) < tigVers) {
      nCompress++;
    }
  }

  if (nErrors > 0) {
    fprintf(stderr, "Store can't be compressed; probably trying to compress to something that isn't the latest version.\n");
    fprintf(stderr, "  "F_U32" tigs failed; "F_U32" compressable\n", nErrors, nCompress);
    delete tigStore;
    exit(1);
  }

  //  Pass 2:  Actually do the moves

  if (nCompress > 0) {
    delete tigStore;
    tigStore = new tgStore(tigName, tigVers, tgStoreModify);
  }

  if (nCompress > 0) {
    fprintf(stderr, "Compressing "F_U32" tigs into version %d\n", nCompress, tigVers);

    for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
      if ((ti % 1000000) == 0)
        fprintf(stderr, "tig %d\n", ti);

      if (tigStore->isDeleted(ti)) {
        continue;
      }

      if (tigStore->getVersion(ti) == tigVers)
        continue;

      tgTig *tig = tigStore->loadTig(ti);

      if (tig == NULL)
        continue;

      tigStore->insertTig(tig, true);
      tigStore->unloadTig(ti);
    }
  }

  //  Now clean up the older files

  if (nCompress > 0) {
    for (uint32 version=1; version<tigVers; version++) {
      fprintf(stderr, "Purge version "F_U32".\n", version);
      tigStore->purgeVersion(version);
    }
  }

  //  And the newer files

  delete tigStore;
}


void
dumpFmap(FILE   *out,
         tgTig  *tig) {
  uint32   fiMax = tig->numberOfChildren();

  for (uint32 fi=0; fi<fiMax; fi++) {
    tgPosition *imp = tig->getChild(fi);

    fprintf(stdout, F_U32"\t"F_U32"\t"F_S32"\t"F_S32"\n",
            imp->ident(), tig->tigID(), imp->bgn(), imp->end());
  }
}




int
main (int argc, char **argv) {
  char          tmpName[FILENAME_MAX] = {0};
  char         *gkpName        = NULL;
  char         *tigName        = NULL;
  int           tigVers        = -1;
  bool          tigIDset       = false;
  uint32        tigIDbgn       = 0;
  uint32        tigIDend       = UINT32_MAX;
  uint32        opType         = 0;
  uint32        dumpFlags      = 0;
  uint32        dumpAll        = 0;
  char         *editName       = NULL;
  char         *replaceName    = NULL;
  bool          sameVersion    = true;
  bool          append         = false;
  char         *buildName      = NULL;

  uint32        minNreads      = 0;
  uint32        maxNreads      = UINT32_MAX;

  uint32        minCoverage    = 0;

  tgTig        *tig             = NULL;

  bool          withQV          = false;
  bool          withDots        = true;
  uint32        displayWidth    = 100;
  uint32        displaySpacing  = 3;

  tgTigSizeAnalysis *siz       = NULL;
  uint64             sizSize   = 0;

  uint64            *cov       = NULL;
  uint64             covMax    = 0;

  char              *outPrefix = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-u") == 0) {
      AS_UTL_decodeRange(argv[++arg], tigIDbgn, tigIDend);
      tigIDset    = true;

    } else if (strcmp(argv[arg], "-U") == 0) {
      dumpAll     = TRUE;
      tigIDset    = true;

    } else if (strcmp(argv[arg], "-d") == 0) {
      arg++;

      opType = OPERATION_TIG;

      if      (strcmp(argv[arg], "properties") == 0)
        dumpFlags = DUMP_PROPERTIES;

      else if (strcmp(argv[arg], "consensus") == 0) {
        dumpFlags = DUMP_CONSENSUS;

        if ((arg+1 < argc) && (isdigit(argv[arg+1][0])))
          minCoverage = atoi(argv[++arg]);
      }

      else if (strcmp(argv[arg], "consensusgapped") == 0)
        dumpFlags = DUMP_CONSENSUSGAPPED;

      else if (strcmp(argv[arg], "layout") == 0)
        dumpFlags = DUMP_LAYOUT;

      else if (strcmp(argv[arg], "multialign") == 0)
        dumpFlags = DUMP_MULTIALIGN;

      else if (strcmp(argv[arg], "sizes") == 0)
        dumpFlags = DUMP_SIZES;

      else if (strcmp(argv[arg], "coverage") == 0)
        dumpFlags = DUMP_COVERAGE;

      else if (strcmp(argv[arg], "overlap") == 0)
        dumpFlags = DUMP_THINOVERLAP;

      else if (strcmp(argv[arg], "fmap") == 0)
        dumpFlags = DUMP_FMAP;

      else
        fprintf(stderr, "%s: Unknown dump option '-d %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      arg++;

      if      (strcmp(argv[arg], "list") == 0)
        opType = OPERATION_LIST;

      else if (strcmp(argv[arg], "properties") == 0)
        opType = OPERATION_PROPERTIES;

      else
        fprintf(stderr, "%s: Unknown dump option '-D %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      opType = OPERATION_EDIT;
      editName = argv[++arg];

    } else if (strcmp(argv[arg], "-A") == 0) {
       append = true;

    } else if (strcmp(argv[arg], "-compress") == 0) {
      opType = OPERATION_COMPRESS;

    } else if (strcmp(argv[arg], "-nreads") == 0) {
      minNreads = atoi(argv[++arg]);
      maxNreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-w") == 0) {
      displayWidth = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      displaySpacing  = atoi(argv[++arg]);
      sizSize         = atoll(argv[arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -G <gkpStore> -T <tigStore> <v> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <gkpStore>         Path to the gatekeeper store\n");
    fprintf(stderr, "  -T <tigStore> <v>     Path to the tigStore, version, to use\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <operation>        Dump something about the store\n");
    fprintf(stderr, "     list               ...a list of the unitigs in the store (NOT IMPLEMENTED)\n");
    fprintf(stderr, "     properties         ...a list of properties for ALL multialigns in the store (for -E)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -u id[-id]            Unitig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive\n");
    fprintf(stderr, "  -c id[-id]            Contig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U                    Dump ALL unitigs (for -d option)\n");
    fprintf(stderr, "  -C                    Dump ALL contigs (for -d option)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nreads min max       Dump tigs with between min and max reads (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d <operation>        Dump something about a multialign (-c or -u) in the store\n");
    fprintf(stderr, "     properties         ...properties\n");
    fprintf(stderr, "     frags              ...a list of fragments\n");
    fprintf(stderr, "     unitigs            ...a list of unitigs\n");
    fprintf(stderr, "     consensus [C]      ...the consensus sequence\n");
    fprintf(stderr, "                             if C supplied, only consensus with coverage >= C is output\n");
    fprintf(stderr, "     consensusgapped    ...the consensus sequence, with gaps as indicated in the multialignment\n");
    fprintf(stderr, "     layout             ...the layout\n");
    fprintf(stderr, "     multialign         ...the full multialignment\n");
    fprintf(stderr, "     sizes              ...an analysis of sizes of the tigs\n");
    fprintf(stderr, "     coverage           ...an analysis of read coverage of the tigs\n");
    fprintf(stderr, "     overlap            ...an analysis of read overlaps in the tigs\n");
    fprintf(stderr, "     fmap               ...a map from fragment IID to unitig IID\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -E <editFile>         Change properties of multialigns\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -compress             Move tigs from earlier versions into the specified version.  This removes\n");
    fprintf(stderr, "                        historical versions of unitigs/contigs, and can save tremendous storage space,\n");
    fprintf(stderr, "                        but makes it impossible to back up the assembly past the specified versions\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d multialign':\n");
    fprintf(stderr, "  -w width              Width of the page.\n");
    fprintf(stderr, "  -s spacing            Spacing between reads on the same line.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d coverage':\n");
    fprintf(stderr, "  -o prefix             Output files will be written to 'prefix.*' in the current directory.\n");
    fprintf(stderr, "                        (defaults to 'tigStore' (the -t option) if not set.)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d sizes':\n");
    fprintf(stderr, "  -s genomesize         Denominator to use for n50 computation\n");
    exit(1);
  }



  if (opType == OPERATION_COMPRESS) {
    operationCompress(tigName, tigVers);
    exit(0);
  }


  gkStore *gkpStore = new gkStore(gkpName);
  tgStore *tigStore = new tgStore(tigName, tigVers);

  if (outPrefix == NULL)
    outPrefix = tigName;


  if ((opType == OPERATION_EDIT) && (editName != NULL)) {
    delete tigStore;
    tigStore = new tgStore(tigName, tigVers, tgStoreModify);
    changeProperties(tigStore, editName);
  }




  if (opType == OPERATION_LIST) {
#warning Not dumping the list of tigs.
    fprintf(stderr, "DUMPING LIST NOT IMPLEMENTED\n");
    //tigStore->dumpMultiAlignRTable(false);
  }


  if (opType == OPERATION_PROPERTIES) {
    for (uint32 i=0; i<tigStore->numTigs(); i++) {
      if (tigStore->isDeleted(i) == false) {
        fprintf(stdout, "coverageStat     %8"F_U32P" %f\n", i, tigStore->getCoverageStat(i));
        //fprintf(stdout, "microhetProb     %8"F_U32P" %f\n", i, tigStore->getMicroHetProb(i));
        fprintf(stdout, "suggestRepeat    %8"F_U32P" %c\n", i, tigStore->getSuggestRepeat(i)   ? 'T' : 'F');
        fprintf(stdout, "suggestUnique    %8"F_U32P" %c\n", i, tigStore->getSuggestUnique(i)   ? 'T' : 'F');
        fprintf(stdout, "suggestCircular  %8"F_U32P" %c\n", i, tigStore->getSuggestCircular(i) ? 'T' : 'F');
        fprintf(stdout, "suggestHaploid   %8"F_U32P" %c\n", i, tigStore->getSuggestHaploid(i)  ? 'T' : 'F');
      }
    }
  }


  if (opType == OPERATION_TIG) {
    if (tigIDset == false) {
      fprintf(stderr, "ERROR: No tig range set with -u or -U.\n");
      exit(1);
    }

    uint32   nTigs = tigStore->numTigs();

    if (dumpAll == TRUE) {
      tigIDbgn = 0;
      tigIDend = nTigs - 1;
    }

    if (nTigs <= tigIDbgn) {
      fprintf(stderr, "ERROR: only "F_U32" tigs in the store (IDs 0-"F_U32" inclusive); can't dump requested range -u "F_U32"-"F_U32"\n",
              nTigs,
              nTigs-1,
              tigIDbgn, tigIDend);
    }

    if (nTigs <= tigIDend)
      tigIDend = nTigs - 1;

    if (dumpFlags == DUMP_SIZES)
      siz = new tgTigSizeAnalysis(sizSize);

    if (dumpFlags == DUMP_COVERAGE) {
      covMax = 1048576;
      cov    = new uint64 [covMax];

      memset(cov, 0, sizeof(uint64) * covMax);
    }

    for (uint32 ti=tigIDbgn; ti<=tigIDend; ti++) {
      uint32  Nreads = tigStore->getNumChildren(ti);

      if ((Nreads    < minNreads) ||
          (maxNreads < Nreads))
        continue;

      tig = tigStore->loadTig(ti);

      if (tig == NULL)
        continue;

      if (dumpFlags == DUMP_PROPERTIES)
        dumpProperties(tigStore, tig);

      if (dumpFlags == DUMP_CONSENSUS)
        dumpConsensus(tigStore, tig, false, minCoverage);

      if (dumpFlags == DUMP_CONSENSUSGAPPED)
        dumpConsensus(tigStore, tig, true, minCoverage);

      if (dumpFlags == DUMP_LAYOUT)
        tig->dumpLayout(stdout);

      if (dumpFlags == DUMP_MULTIALIGN)
        tig->display(stdout, gkpStore, displayWidth, displaySpacing, withQV, withDots);

      if (dumpFlags == DUMP_SIZES)
        siz->evaluateTig(tig);

      if (dumpFlags == DUMP_COVERAGE)
        dumpCoverage(tigStore, tig, 2, UINT32_MAX, cov, covMax, outPrefix);

      if (dumpFlags == DUMP_THINOVERLAP)
        dumpThinOverlap(tigStore, tig, sizSize);

      if (dumpFlags == DUMP_FMAP)
        dumpFmap(stdout, tig);

      tigStore->unloadTig(ti);
    }
  }

  if (siz) {
    siz->finalize();
    siz->printSummary(stdout);
    delete siz;
  }

  if (cov) {
    char  N[FILENAME_MAX];
    FILE *F;

    sprintf(N, "%s.depthHistogram", outPrefix);

    errno = 0;
    F = fopen(N, "w");

    uint32   hMax = covMax - 1;
    while ((hMax > 0) && (cov[hMax] == 0))
      hMax--;

    for (uint32 i=0; i<=hMax; i++)
      fprintf(F, F_U32"\t"F_U64"\n", i, cov[i]);

    fclose(F);

    delete [] cov;

    FILE *gnuPlot = popen("gnuplot > /dev/null 2>&1", "w");

    if (gnuPlot) {
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "set output '%s.depthHistogram.png'\n", outPrefix);
      fprintf(gnuPlot, "set xlabel 'depth'\n");
      fprintf(gnuPlot, "set ylabel 'number of bases'\n");
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "plot '%s.depthHistogram' using 1:2 with lines title '%s tigs %u-%u depthHistogram', \\\n",
              outPrefix, outPrefix, tigIDbgn, tigIDend);

      fclose(gnuPlot);
    }
  }

  delete gkpStore;
  delete tigStore;

  exit(0);
}
