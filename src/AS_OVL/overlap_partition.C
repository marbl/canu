
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

const char *mainid = "$Id: overlap_partition.C,v 1.10 2012-02-21 00:31:58 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_decodeRange.H"
#include "AS_PER_gkpStore.h"


//  Reads gkpStore, outputs four files:
//    ovlbat - batch names
//    ovljob - job names
//    ovlopt - overlapper options

uint32  batchMax = 1000;

void
outputJob(FILE   *BAT,
          FILE   *JOB,
          FILE   *OPT,
          AS_IID  hashBeg,
          AS_IID  hashEnd,
          AS_IID  refBeg,
          AS_IID  refEnd,
          uint32  maxNumFrags,
          uint32  maxLength,
          uint32 &batchSize,
          uint32 &batchName,
          uint32 &jobName) {

  fprintf(BAT, "%03"F_U32P"\n", batchName);
  fprintf(JOB, "%06"F_U32P"\n", jobName);

  if (maxNumFrags == 0) {
    fprintf(OPT, "-h "F_U32"-"F_U32" -r "F_U32"-"F_U32"\n",
            hashBeg, hashEnd, refBeg, refEnd);
    fprintf(stderr, "HASH %10d-%10d  REFR %10d-%10d JOB %d\n",
            hashBeg, hashEnd, refBeg, refEnd, jobName);
  } else {
    fprintf(OPT, "-h "F_U32"-"F_U32" -r "F_U32"-"F_U32" --hashstrings "F_U32" --hashdatalen "F_U32"\n",
            hashBeg, hashEnd, refBeg, refEnd, maxNumFrags, maxLength);
    fprintf(stderr, "HASH %10d-%10d  REFR %10d-%10d  STRINGS %10d  BASES %10d JOB %d\n",
            hashBeg, hashEnd, refBeg, refEnd, maxNumFrags, maxLength, jobName);
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
loadFragmentLengths(gkStore *gkp,
                    set<uint32> &libToHash, AS_IID &hashMin, AS_IID &hashMax,
                    set<uint32> &libToRef,  AS_IID &refMin,  AS_IID &refMax) {
  uint32     numFrags = gkp->gkStore_getNumFragments();
  uint32     numLibs  = gkp->gkStore_getNumLibraries();
  uint32    *fragLen  = new uint32 [numFrags + 1];

  bool testHash = false;
  bool testRef  = false;

  if (libToHash.size() > 0) {
    testHash  = true;
    hashMin   = AS_IID_MAX;
    hashMax   = 0;
  }

  if (libToRef.size() > 0) {
    testRef  = true;
    refMin   = AS_IID_MAX;
    refMax   = 0;
  }

  bool  *doHash = new bool [numLibs + 1];
  bool  *doRef  = new bool [numLibs + 1];

  for (uint32 i=0; i<=numLibs; i++) {
    doHash[i] = (libToHash.count(i) == 0) ? false : true;
    doRef[i]  = (libToRef.count(i)  == 0) ? false : true;
  }

  fprintf(stderr, "Loading lengths of "F_U32" fragments ("F_SIZE_T"mb)\n",
          numFrags, (numFrags * sizeof(uint32)) >> 20);

  memset(fragLen, 0, sizeof(uint32) * (numFrags + 1));

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  for (uint32 ii=1; ii<=numFrags; ii++) {
    fs->next(&fr);

    assert(fr.gkFragment_getReadIID() == ii);

    if (fr.gkFragment_getIsDeleted() == false)
      fragLen[ii] = fr.gkFragment_getClearRegionLength();

    if ((testHash == true) && (doHash[fr.gkFragment_getLibraryIID()] == true)) {
      if (ii < hashMin)
        hashMin = ii;
      if (hashMax < ii)
        hashMax = ii;      
    }

    if ((testRef == true) && (doRef[fr.gkFragment_getLibraryIID()] == true)) {
      if (ii < refMin)
        refMin = ii;
      if (refMax < ii)
        refMax = ii;      
    }

    if ((ii % 1048576) == 0)
      fprintf(stderr, "Loading lengths at "F_U32" out of "F_U32".  H: "F_IID","F_IID"  R: "F_IID","F_IID"\n",
              ii, numFrags, hashMin, hashMax, refMin, refMax);
  }

  delete fs;

  return(fragLen);
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
  AS_IID  hashMin = 1;
  AS_IID  hashBeg = 1;
  AS_IID  hashEnd = 0;
  AS_IID  hashMax = AS_IID_MAX;

  AS_IID  refMin = 1;
  AS_IID  refBeg = 1;
  AS_IID  refEnd = 0;
  AS_IID  refMax = AS_IID_MAX;

  uint32  batchSize = 0;
  uint32  batchName = 1;
  uint32  jobName   = 1;

  uint32     numFrags = gkp->gkStore_getNumFragments();
  uint32    *fragLen  = NULL;

  if ((ovlRefBlockLength > 0) ||
      (libToHash.size() > 0) ||
      (libToRef.size() > 0))
    fragLen = loadFragmentLengths(gkp, libToHash, hashMin, hashMax, libToRef, refMin, refMax);

  if (hashMax > numFrags)
    hashMax = numFrags;
  if (refMax > numFrags)
    refMax = numFrags;

  fprintf(stderr, "Partitioning for hash: "F_IID"-"F_IID" ref: "F_IID","F_IID"\n",
          hashMin, hashMax, refMin, refMax);

  hashBeg = hashMin;

  while (hashBeg < hashMax) {
    hashEnd = hashBeg + ovlHashBlockSize - 1;

    if (hashEnd > hashMax)
      hashEnd = hashMax;

    refBeg = refMin;
    refEnd = 0;

    while ((refBeg < refMax) &&
           (refBeg < hashEnd)) {
      uint64  refLen  = 0;

      if (ovlRefBlockLength > 0) {
        do {
          refEnd++;
          refLen += fragLen[refEnd];
        } while ((refLen < ovlRefBlockLength) && (refEnd < refMax));

      } else {
        refEnd = refBeg + ovlRefBlockSize - 1;
      }

      if (refEnd > refMax)
        refEnd = refMax;
      if (refEnd > hashEnd)
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
  AS_IID  hashMin = 1;
  AS_IID  hashBeg = 1;
  AS_IID  hashEnd = 0;
  AS_IID  hashMax = AS_IID_MAX;

  AS_IID  refMin = 1;
  AS_IID  refBeg = 1;
  AS_IID  refEnd = 0;
  AS_IID  refMax = AS_IID_MAX;

  uint32  batchSize = 0;
  uint32  batchName = 1;
  uint32  jobName   = 1;

  uint32     numFrags = gkp->gkStore_getNumFragments();
  uint32    *fragLen  = loadFragmentLengths(gkp, libToHash, hashMin, hashMax, libToRef, refMin, refMax);

  if (hashMax > numFrags)
    hashMax = numFrags;
  if (refMax > numFrags)
    refMax = numFrags;

  fprintf(stderr, "Partitioning for hash: "F_IID"-"F_IID" ref: "F_IID","F_IID"\n",
          hashMin, hashMax, refMin, refMax);

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  hashBeg = hashMin;
  hashEnd = hashMin - 1;

  while (hashBeg < hashMax) {
    uint64  hashLen = 0;

    assert(hashEnd == hashBeg - 1);

    //  Non deleted fragments contribute one byte per untrimmed base,
    //  and every fragment contributes one more byte for the terminating zero.
    do {
      hashEnd++;
      hashLen += fragLen[hashEnd] + 1;
    } while ((hashLen < ovlHashBlockLength) && (hashEnd < hashMax));

    assert(hashEnd <= hashMax);

    refBeg = refMin;
    refEnd = 0;

    while ((refBeg < refMax) &&
           (refBeg < hashEnd)) {
      uint64  refLen = 0;

      if (ovlRefBlockLength > 0) {
        do {
          refEnd++;
          refLen += fragLen[refEnd];
        } while ((refLen < ovlRefBlockLength) && (refEnd < refMax));

      } else {
        refEnd = refBeg + ovlRefBlockSize - 1;
      }

      if (refEnd > refMax)
        refEnd = refMax;
      if (refEnd > hashEnd)
        refEnd = hashEnd;

      outputJob(BAT, JOB, OPT, hashBeg, hashEnd, refBeg, refEnd, hashEnd - hashBeg + 1, hashLen, batchSize, batchName, jobName);

      refBeg = refEnd + 1;
    }

    hashBeg = hashEnd + 1;
  }

  delete fs;
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

  gkStore   *gkp      = new gkStore(gkpStoreName, FALSE, FALSE, true);

  gkp->gkStore_metadataCaching(true);

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
