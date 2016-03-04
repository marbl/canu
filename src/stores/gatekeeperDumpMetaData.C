
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
#include "gkStore.H"



void
dumpLibs(gkStore *gkp, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping libraries from %u to %u (inclusive).\n", bgnID, endID);

  fprintf(stdout, "libID\tnonRandom\treadType\tcorrectBases\tfinalTrim\tremoveDupe\tremoveSpur\tremoveChimer\tcheckSubRead\tdefaultQV\tlibName\n");

  for (uint32 lid=bgnID; lid<=endID; lid++) {
    gkLibrary  *library = gkp->gkStore_getLibrary(lid);

    fprintf(stdout, F_U32"\t"F_U32"\t%s\t%s\t%s\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t%s\n",
            library->gkLibrary_libraryID(),
            library->gkLibrary_isNonRandom(),
            library->gkLibrary_readTypeString(),
            library->gkLibrary_readCorrectionString(),
            library->gkLibrary_finalTrimString(),
            library->gkLibrary_removeDuplicateReads(),
            library->gkLibrary_removeSpurReads(),
            library->gkLibrary_removeChimericReads(),
            library->gkLibrary_checkForSubReads(),
            library->gkLibrary_defaultQV(),
            library->gkLibrary_libraryName());
  }
}


void
dumpReads(gkStore *gkp, uint32 bgnID, uint32 endID, bool fullDump) {
  //fprintf(stderr, "Dumping reads from %u to %u (inclusive).\n", bgnID, endID);

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    gkRead  *read = gkp->gkStore_getReadInPartition(rid);

    if (read == NULL)
      continue;

    if (fullDump == false)
      fprintf(stdout, F_U32"\t"F_U32"\t"F_U32"\n",
              read->gkRead_readID(),
              read->gkRead_libraryID(),
              read->gkRead_sequenceLength());
    else
      fprintf(stdout, F_U32"\t"F_U32"\t"F_U32"\t"F_U64"\t"F_U64"\n",
              read->gkRead_readID(),
              read->gkRead_libraryID(),
              read->gkRead_sequenceLength(),
              read->gkRead_mPtr(),
              read->gkRead_pID());
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

  void     addRead(gkRead *read) {
    uint32 l = read->gkRead_sequenceLength();

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
dumpStats(gkStore *gkp, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping read statistics from %u to %u (inclusive).\n", bgnID, endID);

  readStats  *rs = new readStats [gkp->gkStore_getNumLibraries() + 1];

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    gkRead  *read = gkp->gkStore_getReadInPartition(rid);

    if (read == NULL)
      continue;

    uint32   l    = read->gkRead_libraryID();

    rs[0].addRead(read);
    rs[l].addRead(read);
  }

  //  Stats per library (this mode when -libs -stats is selected?) and global.
  //  Stats include:
  //    number of reads
  //    total bases
  //    min, mean, stddev, max base per read
  //    length histogram plot

  for (uint32 l=0; l<gkp->gkStore_getNumLibraries() + 1; l++) {
    fprintf(stdout, "library "F_U32"  reads "F_U32" bases: total "F_U64" ave "F_U64" min "F_U64" max "F_U64"\n",
            l, rs[l].numberOfReads(), rs[l].numberOfBases(), rs[l].numberOfBases() / rs[l].numberOfReads(), rs[l].minBases(), rs[l].maxBases());
  }
}


int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  uint32           gkpStorePart      = UINT32_MAX;

  bool             wantLibs          = false;
  bool             wantReads         = true;
  bool             wantStats         = false;  //  Useful only for reads

  bool             fullDump          = true;

  uint32           bgnID             = 1;
  uint32           endID             = UINT32_MAX;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        gkpStorePart = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-libs") == 0) {
      wantLibs  = true;
      wantReads = false;
      wantStats = false;

    } else if (strcmp(argv[arg], "-reads") == 0) {
      wantLibs  = false;
      wantReads = true;
      wantStats = false;

    } else if (strcmp(argv[arg], "-full") == 0) {
      fullDump = true;

    } else if (strcmp(argv[arg], "-stats") == 0) {
      wantLibs  = false;
      wantReads = false;
      wantStats = true;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      bgnID  = atoi(argv[++arg]);
      endID  = bgnID;

    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (gkpStoreName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore [p] [...]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G gkpStore [p]  dump reads from 'gkpStore', restricted to\n");
    fprintf(stderr, "                   partition 'p', if supplied.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -libs            dump information about libraries\n");
    fprintf(stderr, "  -reads [-full]   dump information about reads\n");
    fprintf(stderr, "                     (-full also dumps some storage metadata)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -stats           dump summary statistics on reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b id            output starting at read/library 'id'\n");
    fprintf(stderr, "  -e id            output stopping after read/library 'id'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r id            output only the single read 'id'\n");
    fprintf(stderr, "\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-G) supplied.\n");

    exit(1);
  }

  gkStore    *gkpStore  = gkStore::gkStore_open(gkpStoreName, gkStore_readOnly, gkpStorePart);
  uint32      numReads  = gkpStore->gkStore_getNumReads();
  uint32      numLibs   = gkpStore->gkStore_getNumLibraries();


  if (wantLibs) {
    if (bgnID < 1)         bgnID = 1;
    if (numLibs < endID)   endID = numLibs;

  } else {
    if (bgnID < 1)         bgnID = 1;
    if (numReads < endID)  endID = numReads;
  }


  if (endID < bgnID)
    fprintf(stderr, "No objects to dump; reversed ranges make no sense: bgn="F_U32" end="F_U32"??\n", bgnID, endID);


  if (wantLibs)
    dumpLibs(gkpStore, bgnID, endID);

  if (wantReads)
    dumpReads(gkpStore, bgnID, endID, fullDump);

  if (wantStats)
    dumpStats(gkpStore, bgnID, endID);


  gkpStore->gkStore_close();

  exit(0);
}
