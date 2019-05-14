
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


void
dumpLibs(sqStore *seq, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping libraries from %u to %u (inclusive).\n", bgnID, endID);

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
dumpReads(sqStore *seq, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping reads from %u to %u (inclusive).\n", bgnID, endID);

  fprintf(stdout, "    readID  libraryID     seqLen     rawLen     corLen   clearBgn   clearEnd  segm        byte  part      flags\n");
  fprintf(stdout, "---------- ---------- ---------- ---------- ---------- ---------- ---------- ----- ----------- ----- ----------\n");

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    sqRead  *read = seq->sqStore_getRead(rid);

    if ((read == NULL) ||
        (seq->sqStore_readInPartition(rid) == false))
      continue;

    fprintf(stdout, "%10" F_U32P " %10" F_U32P " %10" F_U32P " %10" F_U32P " %10" F_U32P " %10" F_U32P " %10" F_U32P " %5" F_U64P " %11" F_U64P " %4" F_U64P " %7s%c%c%c\n",
            read->sqRead_readID(),
            read->sqRead_libraryID(),
            read->sqRead_sequenceLength(),
            read->sqRead_sequenceLength(sqRead_raw),
            read->sqRead_sequenceLength(sqRead_corrected),
            read->sqRead_clearBgn(),
            read->sqRead_clearEnd(),
            read->sqRead_mSegm(),
            read->sqRead_mByte(),
            read->sqRead_mPart(),
            "",
            read->sqRead_ignore()  ? 'I' : '-',
            read->sqRead_cExists() ? 'C' : '-',
            read->sqRead_tExists() ? 'T' : '-');
  }
}


class readStats {
public:
  readStats() {
    _nBases   = 0;
    _minBases = UINT32_MAX;
    _maxBases = 0;
  };

  ~readStats() {
  };

  void     addRead(sqRead *read) {
    uint32 l = read->sqRead_sequenceLength();

    _readLengths.push_back(l);

    _nBases += l;

    if (l < _minBases)  _minBases = l;
    if (_maxBases < l)  _maxBases = l;
  };

  uint32   numberOfReads(void)  { return(_readLengths.size());  };
  uint64   numberOfBases(void)  { return(_nBases);              };
  uint64   minBases(void)       { return(_minBases);            };
  uint64   maxBases(void)       { return(_maxBases);            };


private:
  vector<uint32>  _readLengths;

  uint64          _nBases;
  uint32          _minBases;
  uint32          _maxBases;
};



void
dumpStats(sqStore *seq, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping read statistics from %u to %u (inclusive).\n", bgnID, endID);

  readStats  *rs = new readStats [seq->sqStore_getNumLibraries() + 1];

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    sqRead  *read = seq->sqStore_getRead(rid);

    if ((read == NULL) ||
        (seq->sqStore_readInPartition(rid) == false))
      continue;

    uint32   l = read->sqRead_libraryID();

    rs[0].addRead(read);
    rs[l].addRead(read);
  }

  //  Stats per library (this mode when -libs -stats is selected?) and global.
  //  Stats include:
  //    number of reads
  //    total bases
  //    min, mean, stddev, max base per read
  //    length histogram plot

  for (uint32 l=0; l<seq->sqStore_getNumLibraries() + 1; l++)
    fprintf(stdout, "library " F_U32 "  reads " F_U32 " bases: total " F_U64 " ave " F_U64 " min " F_U64 " max " F_U64 "\n",
            l, rs[l].numberOfReads(), rs[l].numberOfBases(), rs[l].numberOfBases() / rs[l].numberOfReads(), rs[l].minBases(), rs[l].maxBases());

  delete [] rs;
}


int
main(int argc, char **argv) {
  char            *seqStoreName      = NULL;
  uint32           seqStorePart      = UINT32_MAX;

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

      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        seqStorePart = atoi(argv[++arg]);

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
    fprintf(stderr, "  -S seqStore [p]  dump reads from 'seqStore', restricted to\n");
    fprintf(stderr, "                   partition 'p', if supplied.\n");
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

  sqStore    *seqStore  = sqStore::sqStore_open(seqStoreName, sqStore_readOnly, seqStorePart);
  uint32      numReads  = seqStore->sqStore_getNumReads();
  uint32      numLibs   = seqStore->sqStore_getNumLibraries();


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


  seqStore->sqStore_close();

  exit(0);
}
