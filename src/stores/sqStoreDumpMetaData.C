
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

#include "runtime.H"
#include "sqStore.H"
#include "strings.H"

#include <map>
#include <vector>


enum class dumpType {
  wantLibs      = 0,
  wantReads     = 1,
  wantStats     = 2,
  wantHistogram = 3,
  wantLengths   = 4,
  wantVersions  = 5,
  wantRevert    = 6
};


void
dumpLibs(sqStore *seq, uint32 bgnID, uint32 endID) {
  uint32  lastLib = seq->sqStore_lastLibraryID();

  assert(bgnID >= 1);
  assert(bgnID <= lastLib);

  assert(endID >= 1);
  assert(endID <= lastLib);

  fprintf(stdout, "ID         technology  name\n");
  fprintf(stdout, "---- ----------------  ------------------------------\n");

  for (uint32 lid=bgnID; lid<=endID; lid++) {
    sqLibrary  *library = seq->sqStore_getLibrary(lid);

    fprintf(stdout, "%-4" F_U32P " %16s  %s\n",
            library->sqLibrary_libraryID(),
            library->sqLibrary_techTypeString(),
            library->sqLibrary_libraryName());
  }
}


void
dumpReads_printHeader(sqRead_which w) {
  char  h1[1024] = {0};
  char  h2[1024] = {0};
  char  h3[1024] = {0};

  if ((w == sqRead_unset) || (((w & sqRead_raw) != 0) && ((w & sqRead_compressed) == 0))) {
    strcat(h1, " --------NORMAL RAW READS--------");
    strcat(h2, "     seqLen   clearBgn   clearEnd");
    strcat(h3, " ---------- ---------- ----------");
  }

  if ((w == sqRead_unset) || (((w & sqRead_raw) != 0) && ((w & sqRead_compressed) != 0))) {
    strcat(h1, " ------COMPRESSED RAW READS------");
    strcat(h2, "     seqLen   clearBgn   clearEnd");
    strcat(h3, " ---------- ---------- ----------");
  }

  if ((w == sqRead_unset) || (((w & sqRead_corrected) != 0) && ((w & sqRead_compressed) == 0))) {
    strcat(h1, " -----NORMAL CORRECTED READS-----");
    strcat(h2, "     seqLen   clearBgn   clearEnd");
    strcat(h3, " ---------- ---------- ----------");
  }

  if ((w == sqRead_unset) || (((w & sqRead_corrected) != 0) && ((w & sqRead_compressed) != 0))) {
    strcat(h1, " ---COMPRESSED CORRECTED READS---");
    strcat(h2, "     seqLen   clearBgn   clearEnd");
    strcat(h3, " ---------- ---------- ----------");
  }

  fprintf(stdout, "                     %s \n", h1);
  fprintf(stdout, "    readID  libraryID%s     flags blob     position\n", h2);
  fprintf(stdout, "---------- ----------%s  -------- ---- ------------\n", h3);
}



bool
dumpReads_setClearString(sqStore *seq, uint32 rid, char *len, char *bgn, char *end, sqRead_which w) {

  //  Set the length of the read.

  if (seq->sqStore_isValidRead(rid, w) == false)
    memcpy(len, "          -", sizeof(char) * 11);

  else if (seq->sqStore_isIgnoredRead(rid, w) == true)
    memcpy(len, "    ignored", sizeof(char) * 11);

  else
    snprintf(len, 12, " %10" F_U32P, seq->sqStore_getReadLength(rid, w));

  assert((w & sqRead_trimmed) == sqRead_unset);   //  Otherwise, length above is trimmed length!

  //  Set the clear range.

  if (seq->sqStore_isTrimmedRead(rid, w) == true) {
    snprintf(bgn, 12, " %10" F_U32P, seq->sqStore_getClearBgn(rid, w | sqRead_trimmed));
    snprintf(end, 12, " %10" F_U32P, seq->sqStore_getClearEnd(rid, w | sqRead_trimmed));
  } else {
    memcpy(bgn, "          -", sizeof(char) * 11);
    memcpy(end, "          -", sizeof(char) * 11);
  }

  len[11] = 0;
  bgn[11] = 0;
  end[11] = 0;

  return(seq->sqStore_isValidRead(rid, w));
}



void
dumpReads_setFlagsString(sqStore *seq, uint32 rid, char *flags) {
  bool   rv = seq->sqStore_isValidRead  (rid, sqRead_raw);
  bool   rt = seq->sqStore_isTrimmedRead(rid, sqRead_raw);
  bool   cv = seq->sqStore_isValidRead  (rid, sqRead_corrected);
  bool   ct = seq->sqStore_isTrimmedRead(rid, sqRead_corrected);

  flags[0] = 'r';   //  Default to non-valid raw and corrected reads.
  flags[1] = '-';
  flags[2] = '-';
  flags[3] = '-';

  flags[4] = 'c';
  flags[5] = '-';
  flags[6] = '-';
  flags[7] = '-';

  if (rv) {
    flags[0] = 'R';
    flags[1] = (seq->sqStore_isIgnoredRead(rid, sqRead_raw)) ? 'I' : 'R';

    if (rt) {
      flags[2] = 'T';
      flags[3] = (seq->sqStore_isIgnoredRead(rid, sqRead_raw | sqRead_trimmed)) ? 'I' : 'T';
    }
  }

  if (cv) {
    flags[4] = 'C';
    flags[5] = (seq->sqStore_isIgnoredRead(rid, sqRead_corrected)) ? 'I' : 'C';

    if (ct) {
      flags[6] = 'T';
      flags[7] = (seq->sqStore_isIgnoredRead(rid, sqRead_corrected | sqRead_trimmed)) ? 'I' : 'T';
    }
  }

  flags[8] = 0;
}



void
dumpReads(sqStore *seq, uint32 bgnID, uint32 endID, sqRead_which w, bool showAll) {
  uint32  lastRead = seq->sqStore_lastReadID();

  assert(bgnID >= 1);
  assert(bgnID <= lastRead);

  assert(endID >= 1);
  assert(endID <= lastRead);

  char   l1[1024] = {0};

  char   s1len[16] = {0}, s1bgn[16] = {0}, s1end[16] = {0};
  char   s2len[16] = {0}, s2bgn[16] = {0}, s2end[16] = {0};
  char   s3len[16] = {0}, s3bgn[16] = {0}, s3end[16] = {0};
  char   s4len[16] = {0}, s4bgn[16] = {0}, s4end[16] = {0};
  char   flags[16] = {0};

  dumpReads_printHeader(w);

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    bool  active = showAll;

    dumpReads_setClearString(seq, rid, s1len, s1bgn, s1end, sqRead_raw);
    dumpReads_setClearString(seq, rid, s2len, s2bgn, s2end, sqRead_raw       | sqRead_compressed);
    dumpReads_setClearString(seq, rid, s3len, s3bgn, s3end, sqRead_corrected);
    dumpReads_setClearString(seq, rid, s4len, s4bgn, s4end, sqRead_corrected | sqRead_compressed);
    dumpReads_setFlagsString(seq, rid, flags);

    l1[0] = 0;

    if ((w == sqRead_unset) || (((w & sqRead_raw) != 0) && ((w & sqRead_compressed) == 0))) {
      active |= seq->sqStore_isValidRead(rid, sqRead_raw);
      strcat(l1, s1len);
      strcat(l1, s1bgn);
      strcat(l1, s1end);
    }

    if ((w == sqRead_unset) || (((w & sqRead_raw) != 0) && ((w & sqRead_compressed) != 0))) {
      active |= seq->sqStore_isValidRead(rid, sqRead_raw);
      strcat(l1, s2len);
      strcat(l1, s2bgn);
      strcat(l1, s2end);
    }

    if ((w == sqRead_unset) || (((w & sqRead_corrected) != 0) && ((w & sqRead_compressed) == 0))) {
      active |= seq->sqStore_isValidRead(rid, sqRead_corrected);
      strcat(l1, s3len);
      strcat(l1, s3bgn);
      strcat(l1, s3end);
    }

    if ((w == sqRead_unset) || (((w & sqRead_corrected) != 0) && ((w & sqRead_compressed) != 0))) {
      active |= seq->sqStore_isValidRead(rid, sqRead_corrected);
      strcat(l1, s4len);
      strcat(l1, s4bgn);
      strcat(l1, s4end);
    }

    if (active)
      fprintf(stdout, "%10" F_U32P " %10" F_U32P "%s  %s %4lu %12lu\n",
              rid,
              seq->sqStore_getLibraryIDForRead(rid),
              l1,
              flags,
              seq->sqStore_getReadSegm(rid),
              seq->sqStore_getReadByte(rid));
  }
}


void
doSummarize_lengthHistogram(uint64               *shortLengths,
                            uint32                shortLengthsLen,
                            std::vector<uint64>  &longLengths,
                            uint64                genomeSize,
                            bool                  limitTo1x);


void
dumpHistogram(sqStore *seq, uint32 bgnID, uint32 endID, bool dumpLengths) {
  std::vector<uint64>   longLengths;
  uint64                nBases = 0;

  //  Build a vector of sequence lengths, pass that to 'sequence' to
  //  generate a pretty picture.

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    uint32  len = seq->sqStore_getReadLength(rid);

    if  (seq->sqStore_isValidRead(rid) == false)        //  Skip invalid reads.
      continue;

    if  (seq->sqStore_isIgnoredRead(rid) == true)       //  Skip ignored reads.
      continue;

    if ((seq->sqStore_isTrimmedRead(rid) == false) &&   //  Skip untrimmed reads,
        (sqRead_defaultVersion & sqRead_trimmed))        //  if we want trimmed reads.
        continue;

    longLengths.push_back(len);

    nBases += len;
  }

  if (dumpLengths == false) {
    char   msg[1024] = {0};

    strcat(msg, "Histogram of");

    if (sqRead_defaultVersion & sqRead_raw)         strcat(msg, " raw");
    if (sqRead_defaultVersion & sqRead_corrected)   strcat(msg, " corrected");

    fprintf(stdout, "Histogram of %s reads:\n", sqRead_getDefaultVersion());

    doSummarize_lengthHistogram(nullptr, 0, longLengths, nBases, false);
  }

  else {
    std::sort(longLengths.begin(), longLengths.end(), std::less<uint64>());

    for (uint64 ii=0; ii<longLengths.size(); ii++)
      fprintf(stdout, "%lu\n", longLengths[ii]);
  }
}



sqStoreInfo
getStats(sqStore *seq, uint32 bgnID, uint32 endID) {
  sqStoreInfo    info;
  sqRead_which   w1 = sqRead_raw;
  sqRead_which   w2 = sqRead_raw       | sqRead_compressed;
  sqRead_which   w3 = sqRead_corrected;
  sqRead_which   w4 = sqRead_corrected | sqRead_compressed;

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    bool  exists = false;

    exists |= info.examineRead(rid, seq->sqStore_getReadSeq(rid, w1), w1);
    exists |= info.examineRead(rid, seq->sqStore_getReadSeq(rid, w2), w2);
    exists |= info.examineRead(rid, seq->sqStore_getReadSeq(rid, w3), w3);
    exists |= info.examineRead(rid, seq->sqStore_getReadSeq(rid, w4), w4);

    if (exists)
      info.sqInfo_addRead();
  }

  return(info);
}



void
dumpStats(sqStore *seq, uint32 bgnID, uint32 endID) {
  getStats(seq, bgnID, endID).writeInfoAsText(stdout);
}



void
dumpVersions(char const *seqStoreName) {
  sqRead_which   w1 = sqRead_raw;
  sqRead_which   w2 = sqRead_raw       | sqRead_trimmed;
  sqRead_which   w3 = sqRead_corrected;
  sqRead_which   w4 = sqRead_corrected | sqRead_trimmed;

  fprintf(stdout, "        |----------- Raw ----------|------- Trimmed Raw ------|-------- Corrected -------|---- Trimmed Corrected ---|\n");
  fprintf(stdout, "Version |     reads          bases |     reads          bases |     reads          bases |     reads          bases |\n");
  fprintf(stdout, "--------|---------- ---------------|---------- ---------------|---------- ---------------|---------- ---------------|\n");

  //  Find the max version available.  Then do some loop weirdness to iterate
  //  through versions 1, 2, 3, 4, ..., 0 (the latest).

  uint32  maxV = sqStore::sqStore_lastVersion(seqStoreName);

  for (uint32 v=1; v <= maxV + 1; v++) {
    uint32       version = (v <= maxV) ? v : 0;

    sqStore     *seqStore = new sqStore(seqStoreName, sqStore_readOnly, version);
    uint32       numReads = seqStore->sqStore_lastReadID();
    uint32       numLibs  = seqStore->sqStore_lastLibraryID();
    sqStoreInfo  info     = getStats(seqStore, 1, numReads);

    delete seqStore;

    fprintf(stdout, "%7s | %9u %14lu | %9u %14lu | %9u %14lu | %9u %14lu |\n",
            (version == 0) ? "latest" : toDec(version),
            info.sqInfo_numReads(w1), info.sqInfo_numBases(w1),
            info.sqInfo_numReads(w2), info.sqInfo_numBases(w2),
            info.sqInfo_numReads(w3), info.sqInfo_numBases(w3),
            info.sqInfo_numReads(w4), info.sqInfo_numBases(w4));
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "If this store is part of a Canu assembly:\n");
  fprintf(stdout, "  Version 1 contains all the input data.\n");
  fprintf(stdout, "  Version 2 contains only maxInputCoverage of the input data.\n");
  fprintf(stdout, "  Version 3 contains corrected reads.\n");
  fprintf(stdout, "  Version 4 contains trimmed reads.\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Versions 3 and 4 will not exist if correction and/or trimming\n");
  fprintf(stdout, "are not computed.\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "The 'latest' version is the current active data.  Numeric versions\n");
  fprintf(stdout, "are historical versions and are not active.\n");
  fprintf(stdout, "\n");
}



int
main(int argc, char **argv) {
  char            *seqStoreName      = NULL;

  sqRead_which     which             = sqRead_unset;
  dumpType         reqDump           = dumpType::wantReads;
  uint32           revertVersion     = 0;
  bool             showAll           = false;

  uint32           bgnID             = 1;
  uint32           endID             = UINT32_MAX;

  uint64           genomeSize        = 0;

  argc = AS_configure(argc, argv, 1);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-libs")      == 0)   { reqDump = dumpType::wantLibs;      }
    else if (strcmp(argv[arg], "-reads")     == 0)   { reqDump = dumpType::wantReads;     }
    else if (strcmp(argv[arg], "-stats")     == 0)   { reqDump = dumpType::wantStats;     }
    else if (strcmp(argv[arg], "-histogram") == 0)   { reqDump = dumpType::wantHistogram; }
    else if (strcmp(argv[arg], "-lengths")   == 0)   { reqDump = dumpType::wantLengths;   }
    else if (strcmp(argv[arg], "-versions")  == 0)   { reqDump = dumpType::wantVersions;  }
    else if (strcmp(argv[arg], "-revert")    == 0)   { reqDump = dumpType::wantRevert;   revertVersion = strtouint32(argv[++arg]); }

    else if (strcmp(argv[arg], "-all") == 0) {
      showAll = true;
    }

    else if (strcmp(argv[arg], "-raw") == 0) {
      which |=  sqRead_raw;
      which &= ~sqRead_corrected;
    }

    else if (strcmp(argv[arg], "-corrected") == 0) {
      which &= ~sqRead_raw;
      which |=  sqRead_corrected;
    }

    else if (strcmp(argv[arg], "-trimmed") == 0) {
      which |=  sqRead_trimmed;
    }

    else if (strcmp(argv[arg], "-untrimmed") == 0) {
      which &= ~sqRead_trimmed;
    }

    else if (strcmp(argv[arg], "-compressed") == 0) {
      which &= ~sqRead_normal;
      which |=  sqRead_compressed;
    }

    else if (strcmp(argv[arg], "-uncompressed") == 0) {
      which |=  sqRead_normal;
      which &= ~sqRead_compressed;
    }

    else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnID, endID);
    }

    else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (seqStoreName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -S seqStore [p] [...]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS:\n");
    fprintf(stderr, "  -S seqStore      dump reads from 'seqStore'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT FORMAT:\n");
    fprintf(stderr, "  -libs            dump information about libraries\n");
    fprintf(stderr, "  -reads           dump information about reads.  the 'read selection' flags will\n");
    fprintf(stderr, "                   restrict the output to just those categories selected.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                     There are four pairs of flags, one for raw, raw-trimmed,\n");
    fprintf(stderr, "                     corrected and corrected-trimmed.  Each pair tells if\n");
    fprintf(stderr, "                     the sequence is valid and if it is ignored.\n");
    fprintf(stderr, "                       1st letter - valid (uppercase) or invalid (lowercase).\n");
    fprintf(stderr, "                       2nd letter - used  (uppercase) or ignored (lowercase).\n");
    fprintf(stderr, "                       1st pair   - raw sequence.\n");
    fprintf(stderr, "                       2nd pair   - raw sequence, trimmed.\n");
    fprintf(stderr, "                       3rd pair   - corrected sequence.\n");
    fprintf(stderr, "                       4th pair   - corrected sequence, trimmed.\n");
    fprintf(stderr, "                     Example:\n");
    fprintf(stderr, "                       RR--c--- - Raw version exists and is used.  Corrected\n");
    fprintf(stderr, "                                  version doesn't exist.\n");
    fprintf(stderr, "                       RR--CCTt - Both raw and corrected versions exist and are used.\n");
    fprintf(stderr, "                                  Corrected trimmed version exists, but is ignored.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -stats           dump summary statistics on reads.  all read types are summarized.\n");
    fprintf(stderr, "                   the 'read selection' options have no effect.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -histogram       dump a length histogram\n");
    fprintf(stderr, "  -lengths         dump sorted read lengths\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -versions        dump a list of the historical metadata saved\n");
    fprintf(stderr, "  -revert n        revert to version 'n' and DESTROY later versions\n");
    fprintf(stderr, "                    ** INCORRECT USE WILL CAUSE GREAT SUFFERING **\n");
    fprintf(stderr, "                    ** (you might as well just rm -rf your asm) **\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION:\n");
    fprintf(stderr, "  Applies to -reads, -histogram and -lengths.  The default for these\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r bgn[-end]     output reads/libraies from `bgn` to `end`, inclusive\n");
    fprintf(stderr, "  -all             show reads not active with the current selection, for -reads\n");
    fprintf(stderr, "  -raw             restrict to 'raw' reads\n");
    fprintf(stderr, "  -corrected       restrict to 'corrected' reads\n");
    fprintf(stderr, "  -untrimmed       restrict to 'untrimmed' reads\n");
    fprintf(stderr, "  -trimmed         restrict to 'trimmed' reads\n");
    fprintf(stderr, "  -uncompressed    restrict to 'normal non-homopolymer compressed' reads\n");
    fprintf(stderr, "  -compressed      restrict to 'homopolymer compressed' reads\n");
    fprintf(stderr, "\n");

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: no seqStore (-S) supplied.\n");

    exit(1);
  }

  //  Set the default version if one was supplied.  If one isn't supplied,
  //  which is sqRead_unset, and this does nothing.
  //
  if (which != sqRead_unset)
    sqRead_setDefaultVersion(which);

  sqStore    *seqStore  = new sqStore(seqStoreName, sqStore_readOnly);
  uint32      numReads  = seqStore->sqStore_lastReadID();
  uint32      numLibs   = seqStore->sqStore_lastLibraryID();

  //  If there was no default version supplied, set it to the current default
  //  (as set by the store).  This lets us call sqRead_length et al. with
  //  this 'which' regardless of if it's set by the user or not (otherwise
  //  we'd call sqRead_length(sqRead_unset) which is invalid).
  //
  if (which == sqRead_unset)
    which = sqRead_defaultVersion;

  //fprintf(stderr, "Opened seqStore '%s' for '%s' reads.\n", seqStoreName, sqRead_getDefaultVersion());

  if (reqDump == dumpType::wantLibs) {
    if (bgnID < 1)         bgnID = 1;
    if (numLibs < endID)   endID = numLibs;

  } else {
    if (bgnID < 1)         bgnID = 1;
    if (numReads < endID)  endID = numReads;
  }

  if (endID < bgnID) {
    fprintf(stderr, "ERROR: No objects to dump; reversed ranges make no sense: bgn=" F_U32 " end=" F_U32 ".\n", bgnID, endID);
    exit(1);
  }

  if (reqDump == dumpType::wantLibs)
    dumpLibs(seqStore, bgnID, endID);

  if (reqDump == dumpType::wantReads)
    dumpReads(seqStore, bgnID, endID, which, showAll);

  if (reqDump == dumpType::wantStats)
    dumpStats(seqStore, bgnID, endID);

  if (reqDump == dumpType::wantHistogram)
    dumpHistogram(seqStore, bgnID, endID, false);

  if (reqDump == dumpType::wantLengths)
    dumpHistogram(seqStore, bgnID, endID, true);

  if (reqDump == dumpType::wantVersions)
    dumpVersions(seqStoreName);

  if ((reqDump == dumpType::wantRevert) &&
      (revertVersion != 0))
    sqStore::sqStore_revertVersion(seqStoreName, revertVersion);

  delete seqStore;

  exit(0);
}
