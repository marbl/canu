
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
 *    src/stores/gatekeeperDumpMetaData.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-MAR-18 to 2015-SEP-21
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "sqStore.H"
#include "strings.H"

#include <map>

using namespace std;


void
dumpLibs(sqStore *seq, uint32 bgnID, uint32 endID) {

  fprintf(stdout, "libID\tnonRandom\treadType\tcorrectBases\tfinalTrim\tremoveDupe\tremoveSpur\tremoveChimer\tcheckSubRead\tdefaultQV\tlibName\n");

  for (uint32 lid=bgnID; lid<=endID; lid++) {
    sqLibrary  *library = seq->sqStore_getLibrary(lid);

    fprintf(stdout, F_U32"\t" F_U32 "\t%s\t%s\t%s\t" F_U32 "\t" F_U32 "\t" F_U32 "\t" F_U32 "\t" F_U32 "\t%s\n",
            library->sqLibrary_libraryID(),
            library->sqLibrary_isNonRandom(),
            library->sqLibrary_readTypeString(),
            library->sqLibrary_readCorrectionString(),
            library->sqLibrary_finalTrimString(),
            library->sqLibrary_removeDuplicateReads(),
            library->sqLibrary_removeSpurReads(),
            library->sqLibrary_removeChimericReads(),
            library->sqLibrary_checkForSubReads(),
            library->sqLibrary_defaultQV(),
            library->sqLibrary_libraryName());
  }
}


void
dumpReads_setClearString(sqStore *seqs, uint32 rid, char *len, char *bgn, char *end, sqRead_which w) {

  if (seqs->sqStore_isValidRead(rid, w) == false)
    memcpy(len, "         -", sizeof(char) * 10);

  else if (seqs->sqStore_isIgnoredRead(rid, w) == true)
    memcpy(len, "   ignored", sizeof(char) * 10);

  else
    snprintf(len, 11, "%10" F_U32P, seqs->sqStore_getReadLength(rid, w));

  assert((w & sqRead_trimmed) == sqRead_unset);   //  Otherwise, length above is trimmed length!

  if (seqs->sqStore_isTrimmedRead(rid, w) == true) {
    snprintf(bgn, 11, "%10" F_U32P, seqs->sqStore_getClearBgn(rid, w));
    snprintf(end, 11, "%10" F_U32P, seqs->sqStore_getClearEnd(rid, w));
  } else {
    memcpy(bgn, "         -", sizeof(char) * 10);
    memcpy(end, "         -", sizeof(char) * 10);
  }

  len[10] = 0;
  bgn[10] = 0;
  end[10] = 0;
}



void
dumpReads(sqStore *seqs, uint32 bgnID, uint32 endID) {
  char   s1len[16] = {0}, s1bgn[16] = {0}, s1end[16] = {0};
  char   s2len[16] = {0}, s2bgn[16] = {0}, s2end[16] = {0};
  char   s3len[16] = {0}, s3bgn[16] = {0}, s3end[16] = {0};
  char   s4len[16] = {0}, s4bgn[16] = {0}, s4end[16] = {0};

  fprintf(stdout, "                      --------NORMAL RAW READS-------- ------COMPRESSED RAW READS------ -----NORMAL CORRECTED READS----- ---COMPRESSED CORRECTED READS--- \n");
  fprintf(stdout, "    readID  libraryID     seqLen   clearBgn   clearEnd     seqLen   clearBgn   clearEnd     seqLen   clearBgn   clearEnd     seqLen   clearBgn   clearEnd blobFile    blobPos\n");
  fprintf(stdout, "---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- -------- ----------\n");

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    dumpReads_setClearString(seqs, rid, s1len, s1bgn, s1end, sqRead_raw);
    dumpReads_setClearString(seqs, rid, s2len, s2bgn, s2end, sqRead_raw       | sqRead_compressed);
    dumpReads_setClearString(seqs, rid, s3len, s3bgn, s3end, sqRead_corrected);
    dumpReads_setClearString(seqs, rid, s4len, s4bgn, s4end, sqRead_corrected | sqRead_compressed);

    fprintf(stdout, "%10" F_U32P " %10" F_U32P " %s %s %s %s %s %s %s %s %s %s %s %s %8" F_U64P " %10" F_U64P "\n",
            rid,
            seqs->sqStore_getLibraryIDForRead(rid),
            s1len, s1bgn, s1end,
            s2len, s2bgn, s2end,
            s3len, s3bgn, s3end,
            s4len, s4bgn, s4end,
            seqs->sqStore_getMeta(rid)->sqRead_mSegm(),
            seqs->sqStore_getMeta(rid)->sqRead_mByte());
  }
}



void
dumpStats(sqStore *seqs, uint32 bgnID, uint32 endID) {
  sqStoreInfo    info;
  sqRead_which   w1 = sqRead_raw;
  sqRead_which   w2 = sqRead_raw       | sqRead_compressed;
  sqRead_which   w3 = sqRead_corrected;
  sqRead_which   w4 = sqRead_corrected | sqRead_compressed;

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    info.examineRead(seqs->sqStore_getReadSeq(rid, w1), w1);
    info.examineRead(seqs->sqStore_getReadSeq(rid, w2), w2);
    info.examineRead(seqs->sqStore_getReadSeq(rid, w3), w3);
    info.examineRead(seqs->sqStore_getReadSeq(rid, w4), w4);
  }

  info.writeInfoAsText(stdout);
}




int
main(int argc, char **argv) {
  char            *seqStoreName      = NULL;

  bool             wantLibs          = false;
  bool             wantReads         = true;
  bool             wantStats         = false;  //  Useful only for reads

  uint32           bgnID             = 1;
  uint32           endID             = UINT32_MAX;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-libs") == 0) {
      wantLibs  = true;
      wantReads = false;
      wantStats = false;

    } else if (strcmp(argv[arg], "-reads") == 0) {
      wantLibs  = false;
      wantReads = true;
      wantStats = false;

    } else if (strcmp(argv[arg], "-stats") == 0) {
      wantLibs  = false;
      wantReads = false;
      wantStats = true;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnID, endID);

    } else {
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
    fprintf(stderr, "  -S seqStore      dump reads from 'seqStore'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -libs            dump information about libraries\n");
    fprintf(stderr, "  -reads           dump information about reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -stats           dump summary statistics on reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b id            output starting at read/library 'id'     DEPRECATED\n");
    fprintf(stderr, "  -e id            output stopping after read/library 'id'  DEPRECATED");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r bgn[-end]     output reads/libraies from `bgn` to `end`, inclusive\n");
    fprintf(stderr, "\n");

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: no seqStore (-S) supplied.\n");

    exit(1);
  }

  sqStore    *seqStore  = new sqStore(seqStoreName, sqStore_readOnly);
  uint32      numReads  = seqStore->sqStore_lastReadID();
  uint32      numLibs   = seqStore->sqStore_lastLibraryID();


  if (wantLibs) {
    if (bgnID < 1)         bgnID = 1;
    if (numLibs < endID)   endID = numLibs;

  } else {
    if (bgnID < 1)         bgnID = 1;
    if (numReads < endID)  endID = numReads;
  }


  if (endID < bgnID)
    fprintf(stderr, "No objects to dump; reversed ranges make no sense: bgn=" F_U32 " end=" F_U32 "??\n", bgnID, endID);


  if (wantLibs)
    dumpLibs(seqStore, bgnID, endID);

  if (wantReads)
    dumpReads(seqStore, bgnID, endID);

  if (wantStats)
    dumpStats(seqStore, bgnID, endID);


  delete seqStore;

  exit(0);
}
