
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: gatekeeperDumpFASTQ.C 6775 2015-03-18 03:55:15Z bri $";

#include "AS_global.H"
#include "gkStore.H"



void
dumpLibs(gkStore *gkp, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping libraries from %u to %u (inclusive).\n", bgnID, endID);

  for (uint32 lid=bgnID; lid<=endID; lid++) {
    gkLibrary  *library = gkp->gkStore_getLibrary(lid);

    fprintf(stdout, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t%s\n",
            library->gkLibrary_libraryID(),
            library->gkLibrary_isNonRandom(),
            library->gkLibrary_trustHomopolymerRuns(),
            library->gkLibrary_correctBases(),
            library->gkLibrary_finalTrim(),
            library->gkLibrary_removeDuplicateReads(),
            library->gkLibrary_removeSpurReads(),
            library->gkLibrary_removeChimericReads(),
            library->gkLibrary_checkForSubReads(),
            library->gkLibrary_defaultQV(),
            library->gkLibrary_libraryName());
  }
}


void
dumpReads(gkStore *gkp, uint32 bgnID, uint32 endID) {
  //fprintf(stderr, "Dumping reads from %u to %u (inclusive).\n", bgnID, endID);

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    gkRead  *read = gkp->gkStore_getReadInPartition(rid);

    if (read == NULL)
      continue;

    fprintf(stdout, F_U32"\t"F_U32"\t"F_U32"\t%c\t"F_U64"\t"F_U64"\n",
            read->gkRead_readID(),
            read->gkRead_libraryID(),
            read->gkRead_sequenceLength(),
            read->gkRead_isDeleted() ? 'D' : 'v',
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
    fprintf(stdout, "library "F_U32"  reads "F_U32" bases: total "F_U64" ave "F_U32" min "F_U32" max "F_U32"\n",
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

  uint32           bgnID             = 1;
  uint32           endID             = AS_MAX_READS;

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
    fprintf(stderr, "  -reads           dump information about reads\n");
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

  gkStore    *gkpStore  = new gkStore(gkpStoreName, gkStore_readOnly, gkpStorePart);
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
    dumpReads(gkpStore, bgnID, endID);

  if (wantStats)
    dumpStats(gkpStore, bgnID, endID);


  delete gkpStore;

  exit(0);
}
