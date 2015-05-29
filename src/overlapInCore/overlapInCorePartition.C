
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2011, J. Craig Venter Institute.
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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "gkStore.H"
#include "AS_UTL_decodeRange.H"

//  Reads gkpStore, outputs three files:
//    ovlbat - batch names
//    ovljob - job names
//    ovlopt - overlapper options

uint32  batchMax = 1000;

void
outputJob(FILE   *BAT,
          FILE   *JOB,
          FILE   *OPT,
          uint32  hashBeg,
          uint32  hashEnd,
          uint32  refBeg,
          uint32  refEnd,
          uint32  maxNumReads,
          uint32  maxLength,
          uint32 &batchSize,
          uint32 &batchName,
          uint32 &jobName) {

  fprintf(BAT, "%03"F_U32P"\n", batchName);
  fprintf(JOB, "%06"F_U32P"\n", jobName);

  if (maxNumReads == 0) {
    fprintf(OPT, "-h "F_U32"-"F_U32" -r "F_U32"-"F_U32"\n",
            hashBeg, hashEnd, refBeg, refEnd);
    fprintf(stderr, "HASH %10d-%10d  REFR %10d-%10d JOB %d\n",
            hashBeg, hashEnd, refBeg, refEnd, jobName);
  } else {
    fprintf(OPT, "-h "F_U32"-"F_U32" -r "F_U32"-"F_U32" --hashstrings "F_U32" --hashdatalen "F_U32"\n",
            hashBeg, hashEnd, refBeg, refEnd, maxNumReads, maxLength);
    fprintf(stderr, "HASH %10d-%10d  REFR %10d-%10d  STRINGS %10d  BASES %10d JOB %d\n",
            hashBeg, hashEnd, refBeg, refEnd, maxNumReads, maxLength, jobName);
  }
  refBeg = refEnd + 1;

  batchSize++;

  if (batchSize >= batchMax) {
    batchSize = 0;
    batchName++;
  }

  jobName++;
}



uint32 *
loadReadLengths(gkStore *gkp,
                    set<uint32> &libToHash, uint32 &hashMin, uint32 &hashMax,
                    set<uint32> &libToRef,  uint32 &refMin,  uint32 &refMax) {
  uint32     numReads = gkp->gkStore_getNumReads();
  uint32     numLibs  = gkp->gkStore_getNumLibraries();
  uint32    *readLen  = new uint32 [numReads + 1];

  bool testHash = false;
  bool testRef  = false;

  if (libToHash.size() > 0) {
    testHash  = true;
    hashMin   = UINT32_MAX;
    hashMax   = 0;
  }

  if (libToRef.size() > 0) {
    testRef  = true;
    refMin   = UINT32_MAX;
    refMax   = 0;
  }

  bool  *doHash = new bool [numLibs + 1];
  bool  *doRef  = new bool [numLibs + 1];

  for (uint32 i=0; i<=numLibs; i++) {
    doHash[i] = (libToHash.count(i) == 0) ? false : true;
    doRef[i]  = (libToRef.count(i)  == 0) ? false : true;
  }

  fprintf(stderr, "Loading lengths of "F_U32" fragments ("F_SIZE_T"mb)\n",
          numReads, (numReads * sizeof(uint32)) >> 20);

  memset(readLen, 0, sizeof(uint32) * (numReads + 1));

  for (uint32 ii=1; ii<=numReads; ii++) {
    gkRead  *read = gkp->gkStore_getRead(ii);

    if (read->gkRead_readID() != ii)
      fprintf(stderr, "ERROR: readID=%u != ii=%u\n",
              read->gkRead_readID(), ii);
    assert(read->gkRead_readID() == ii);

    readLen[ii] = read->gkRead_sequenceLength();

    if ((testHash == true) && (doHash[read->gkRead_libraryID()] == true)) {
      if (ii < hashMin)
        hashMin = ii;
      if (hashMax < ii)
        hashMax = ii;      
    }

    if ((testRef == true) && (doRef[read->gkRead_libraryID()] == true)) {
      if (ii < refMin)
        refMin = ii;
      if (refMax < ii)
        refMax = ii;      
    }

    if ((ii % 1048576) == 0)
      fprintf(stderr, "Loading lengths at "F_U32" out of "F_U32".  H: "F_U32","F_U32"  R: "F_U32","F_U32"\n",
              ii, numReads, hashMin, hashMax, refMin, refMax);
  }

  return(readLen);
}


void
partitionFrags(gkStore      *gkp,
               FILE         *BAT,
               FILE         *JOB,
               FILE         *OPT,
               uint64        ovlHashBlockSize,
               uint64        ovlRefBlockLength,
               uint64        ovlRefBlockSize,
               set<uint32>  &libToHash,
               set<uint32>  &libToRef) {
  uint32  hashMin = 1;
  uint32  hashBeg = 1;
  uint32  hashEnd = 0;
  uint32  hashMax = UINT32_MAX;

  uint32  refMin = 1;
  uint32  refBeg = 1;
  uint32  refEnd = 0;
  uint32  refMax = UINT32_MAX;

  uint32  batchSize = 0;
  uint32  batchName = 1;
  uint32  jobName   = 1;

  uint32  numReads = gkp->gkStore_getNumReads();
  uint32 *readLen  = NULL;

  if ((ovlRefBlockLength > 0) ||
      (libToHash.size() > 0) ||
      (libToRef.size() > 0))
    readLen = loadReadLengths(gkp, libToHash, hashMin, hashMax, libToRef, refMin, refMax);

  if (hashMax > numReads)
    hashMax = numReads;
  if (refMax > numReads)
    refMax = numReads;

  fprintf(stderr, "Partitioning for hash: "F_U32"-"F_U32" ref: "F_U32","F_U32"\n",
          hashMin, hashMax, refMin, refMax);

  hashBeg = hashMin;

  while (hashBeg < hashMax) {
    hashEnd = hashBeg + ovlHashBlockSize - 1;

    if (hashEnd > hashMax)
      hashEnd = hashMax;

    refBeg = refMin;
    refEnd = 0;

    while ((refBeg < refMax) &&
           ((refBeg < hashEnd) || (libToHash.size() != 0 && libToHash == libToRef))) {
      uint64  refLen  = 0;

      if (ovlRefBlockLength > 0) {
        do {
          refEnd++;
          refLen += readLen[refEnd];
        } while ((refLen < ovlRefBlockLength) && (refEnd < refMax));

      } else {
        refEnd = refBeg + ovlRefBlockSize - 1;
      }

      if (refEnd > refMax)
        refEnd = refMax;
      if ((refEnd > hashEnd) && (libToHash.size() == 0 || libToHash != libToRef))
        refEnd = hashEnd;

      outputJob(BAT, JOB, OPT, hashBeg, hashEnd, refBeg, refEnd, 0, 0, batchSize, batchName, jobName);

      refBeg = refEnd + 1;
    }

    hashBeg = hashEnd + 1;
  }
}





void
partitionLength(gkStore      *gkp,
                FILE         *BAT,
                FILE         *JOB,
                FILE         *OPT,
                uint64        ovlHashBlockLength,
                uint64        ovlRefBlockLength,
                uint64        ovlRefBlockSize,
                set<uint32>  &libToHash,
                set<uint32>  &libToRef) {
  uint32  hashMin = 1;
  uint32  hashBeg = 1;
  uint32  hashEnd = 0;
  uint32  hashMax = UINT32_MAX;

  uint32  refMin = 1;
  uint32  refBeg = 1;
  uint32  refEnd = 0;
  uint32  refMax = UINT32_MAX;

  uint32  batchSize = 0;
  uint32  batchName = 1;
  uint32  jobName   = 1;

  uint32  numReads = gkp->gkStore_getNumReads();
  uint32 *readLen  = loadReadLengths(gkp, libToHash, hashMin, hashMax, libToRef, refMin, refMax);

  if (hashMax > numReads)
    hashMax = numReads;
  if (refMax > numReads)
    refMax = numReads;

  fprintf(stderr, "Partitioning for hash: "F_U32"-"F_U32" ref: "F_U32","F_U32"\n",
          hashMin, hashMax, refMin, refMax);

  hashBeg = hashMin;
  hashEnd = hashMin - 1;

  while (hashBeg < hashMax) {
    uint64  hashLen = 0;

    assert(hashEnd == hashBeg - 1);

    //  Non deleted reads contribute one byte per untrimmed base, and every fragment contributes one
    //  more byte for the terminating zero.  In 3g, there are no deleted reads.

    do {
      hashEnd++;
      hashLen += readLen[hashEnd] + 1;
    } while ((hashLen < ovlHashBlockLength) && (hashEnd < hashMax));

    assert(hashEnd <= hashMax);

    refBeg = refMin;
    refEnd = 0;

    while ((refBeg < refMax) &&
           ((refBeg < hashEnd) || (libToHash.size() != 0 && libToHash == libToRef))) {
      uint64  refLen = 0;

      if (ovlRefBlockLength > 0) {
        do {
          refEnd++;
          refLen += readLen[refEnd];
        } while ((refLen < ovlRefBlockLength) && (refEnd < refMax));

      } else {
        refEnd = refBeg + ovlRefBlockSize - 1;
      }

      if (refEnd > refMax)
        refEnd = refMax;
      if ((refEnd > hashEnd) && (libToHash.size() == 0 || libToHash != libToRef))
        refEnd = hashEnd;

      outputJob(BAT, JOB, OPT, hashBeg, hashEnd, refBeg, refEnd, hashEnd - hashBeg + 1, hashLen, batchSize, batchName, jobName);

      refBeg = refEnd + 1;
    }

    hashBeg = hashEnd + 1;
  }
}




int
main(int argc, char **argv) {
  char            *gkpStoreName        = NULL;
  gkStore         *gkpStore            = NULL;

  char            *outputPrefix        = NULL;
  char             outputName[FILENAME_MAX];

  uint64           ovlHashBlockLength  = 0;
  uint64           ovlHashBlockSize    = 0;
  uint64           ovlRefBlockLength   = 0;
  uint64           ovlRefBlockSize     = 0;
  bool             checkAllLibUsed     = true;

  set<uint32>      libToHash;
  set<uint32>      libToRef;

  AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-bl") == 0) {
      ovlHashBlockLength = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-bs") == 0) {
      ovlHashBlockSize   = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-rl") == 0) {
      ovlRefBlockLength  = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-rs") == 0) {
      ovlRefBlockSize    = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-H") == 0) {
      AS_UTL_decodeRange(argv[++arg], libToHash);

    } else if (strcmp(argv[arg], "-R") == 0) {
      AS_UTL_decodeRange(argv[++arg], libToRef);

    } else if (strcmp(argv[arg], "-C") == 0) {
       checkAllLibUsed = false;

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else {
      fprintf(stderr, "ERROR:  Unknown option '%s'\n", arg[argv]);
      err++;
    }

    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    exit(1);
  }

  if ((ovlHashBlockLength > 0) && (ovlHashBlockSize > 0))
    fprintf(stderr, "ERROR:  At most one of -bl and -bs can be non-zero.\n"), exit(1);

  if ((ovlRefBlockLength > 0) && (ovlRefBlockSize > 0))
    fprintf(stderr, "ERROR:  At most one of -rl and -rs can be non-zero.\n"), exit(1);

  fprintf(stderr, "HASH: "F_U64" reads or "F_U64" length.\n", ovlHashBlockSize, ovlHashBlockLength);
  fprintf(stderr, "REF:  "F_U64" reads or "F_U64" length.\n", ovlRefBlockSize,  ovlRefBlockLength);

  gkStore   *gkp         = new gkStore(gkpStoreName);
  uint32     numLibs     = gkp->gkStore_getNumLibraries();
  uint32     invalidLibs = 0;

  for (set<uint32>::iterator it=libToHash.begin(); it != libToHash.end(); it++)
    if (numLibs < *it)
      fprintf(stderr, "ERROR: -H "F_U32" is invalid; only "F_U32" libraries in '%s'\n",
              *it, numLibs, gkpStoreName), invalidLibs++;

  for (set<uint32>::iterator it=libToRef.begin(); it != libToRef.end(); it++)
    if (numLibs < *it)
      fprintf(stderr, "ERROR: -R "F_U32" is invalid; only "F_U32" libraries in '%s'\n",
              *it, numLibs, gkpStoreName), invalidLibs++;

  if ((libToHash.size() > 0) && (libToRef.size() > 0)) {
    for (uint32 lib=1; lib<=numLibs; lib++) {
      if ((libToHash.find(lib) == libToHash.end()) &&
          (libToRef.find(lib)  == libToRef.end())) {
        if (checkAllLibUsed == true)
          fprintf(stderr, "ERROR: library "F_U32" is not mentioned in either -H or -R.\n", lib), invalidLibs++;
        else
          fprintf(stderr, "Warning: library "F_U32" is not mentioned in either -H or -R.\n", lib);
       }
    }
  }
  if (invalidLibs > 0)
    fprintf(stderr, "ERROR: one of -H and/or -R are invalid.\n"), exit(1);

  errno = 0;

  sprintf(outputName, "%s/ovlbat", outputPrefix);
  FILE *BAT = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s/ovljob", outputPrefix);
  FILE *JOB = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s/ovlopt", outputPrefix);
  FILE *OPT = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outputName, strerror(errno)), exit(1);

  if (ovlHashBlockLength == 0)
    partitionFrags(gkp, BAT, JOB, OPT, ovlHashBlockSize, ovlRefBlockLength, ovlRefBlockSize, libToHash, libToRef);
  else
    partitionLength(gkp, BAT, JOB, OPT, ovlHashBlockLength, ovlRefBlockLength, ovlRefBlockSize, libToHash, libToRef);

  fclose(BAT);
  fclose(JOB);
  fclose(OPT);

  delete gkp;

  exit(0);
}
