
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

#include "AS_PER_gkpStore.H"

#include "overlapStore.H"

#include "AS_OVS_overlap.H"
#include "AS_OVS_overlapFile.H"
#include "AS_OVS_overlapStore.H"
#include "AS_OBT_acceptableOverlap.H"

#include <ctype.h>
#include <unistd.h>  //  sysconf()

#include <vector>
#include <algorithm>

using namespace std;

#define AS_OVS_CURRENT_VERSION  2

#undef  DELETE_INTERMEDIATE_EARLY
#define DELETE_INTERMEDIATE_LATE


//  This should be private to AS_OVS
//
void
writeOverlaps(char                *ovlName,
              OVSoverlap          *overlapsort,
              uint64               numOvl,
              uint32               jobIndex) {

	char                        name[FILENAME_MAX];

  uint32                      currentFileIndex = jobIndex;
  uint64                      overlapsThisFile = 0;

	OverlapStoreOffsetRecord    offset; 
  OverlapStoreOffsetRecord    missing;
 
  offset.a_iid     = overlapsort[0].a_iid;
  offset.numOlaps  = 0;
	offset.fileno    = jobIndex;

	missing.a_iid    = overlapsort[0].a_iid;
  missing.numOlaps = 0;
  missing.fileno   = jobIndex;


	OverlapStoreInfo ovs;

	ovs.ovsMagic              = 1;
	ovs.ovsVersion            = AS_OVS_CURRENT_VERSION;
  ovs.numOverlapsPerFile    = 1024 * 1024 * 1024 / sizeof(OVSoverlapINT);
  ovs.smallestIID           = UINT64_MAX;
  ovs.largestIID            = 0;
  ovs.numOverlapsTotal      = 0;
  ovs.highestFileIndex      = 0;
	ovs.maxReadLenInBits      = AS_READ_MAX_NORMAL_LEN_BITS;

  sprintf(name, "%s/%04d", ovlName, jobIndex);
  BinaryOverlapFile *bof = AS_OVS_createBinaryOverlapFile(name, TRUE);

	sprintf(name,"%s/%04d.idx", ovlName, jobIndex);

  errno = 0;
  FILE *offsetFile=fopen(name,"w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  fprintf(stderr, "Writing "F_U64" overlaps.\n", numOvl);

	for (uint64 i=0; i<numOvl; i++ ) {
    AS_OVS_writeOverlap(bof, overlapsort + i);

    if (offset.a_iid > overlapsort[i].a_iid) {
			fprintf(stderr, "LAST:  a:"F_U32"\n", offset.a_iid);
			fprintf(stderr, "THIS:  a:"F_U32" b:"F_U32"\n", overlapsort[i].a_iid, overlapsort[i].b_iid);
		}
    assert(offset.a_iid <= overlapsort[i].a_iid);

    ovs.smallestIID = MIN(ovs.smallestIID, overlapsort[i].a_iid);
    ovs.largestIID  = MAX(ovs.largestIID,  overlapsort[i].a_iid);

		//  Put the index to disk, filling any gaps
		if ((offset.numOlaps != 0) && (offset.a_iid != overlapsort[i].a_iid)) {
			while (missing.a_iid < offset.a_iid) {
				missing.fileno    = offset.fileno;
				missing.offset    = offset.offset;
				missing.numOlaps  = 0;

				AS_UTL_safeWrite(offsetFile, &missing, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);
				missing.a_iid++;
			}

			//  One more, since this iid is not missing -- we write it next!
			missing.a_iid++;

			AS_UTL_safeWrite(offsetFile, &offset, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);
			offset.numOlaps  = 0;
		}

		//  Update the index if this is the first overlap for this a_iid
		if (offset.numOlaps == 0) {
			offset.a_iid     = overlapsort[i].a_iid;
			offset.fileno    = currentFileIndex;
			offset.offset    = overlapsThisFile;
		}

		offset.numOlaps++;

		ovs.numOverlapsTotal++;

		overlapsThisFile++;
	}

  //write final a_iid index

	while (missing.a_iid < offset.a_iid) {
		missing.fileno    = offset.fileno;
		missing.offset    = offset.offset;
		missing.numOlaps  = 0;

		AS_UTL_safeWrite(offsetFile, &missing, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);
		missing.a_iid++;
	}

	AS_UTL_safeWrite(offsetFile, &offset, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);

	fclose(offsetFile);

  //  In the nasty case that there were no overlaps in this slice, set meaningful smallest and
  //  largest.  Well, at least, set non-nonsense smallest and largest.

  if (overlapsThisFile == 0) {
    ovs.smallestIID = 0;
    ovs.largestIID = 0;
  }

  {
    sprintf(name,"%s/%04d.ovs", ovlName, jobIndex);

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);
  
    AS_UTL_safeWrite(F, &ovs, "Partition ovs file", sizeof(OverlapStoreInfo), 1);

    fclose(F);	

    fprintf(stderr, "Wrote "F_U64" overlaps into '%s'\n", ovs.numOverlapsTotal, name);
    fprintf(stderr, "Smallest "F_U64" largest "F_U64"\n", ovs.smallestIID, ovs.largestIID);
  }

  AS_OVS_closeBinaryOverlapFile(bof);    
}
















int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  uint32          fileLimit    = 512;   //  Number of 'slices' from bucketizer
  uint32          jobIndex     = 0;     //  'slice' that we are going to be sorting
  uint32          jobIdxMax    = 0;     //  Number of 'buckets' from bucketizer

  uint64          maxMemory    = UINT64_MAX;

  bool            deleteIntermediateEarly = false;
  bool            deleteIntermediateLate  = false;

  bool            forceRun = false;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      //unused gkpName = argv[++arg];
      ++arg;

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-job") == 0) {
      jobIndex  = atoi(argv[++arg]);
      jobIdxMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      maxMemory  = atoi(argv[++arg]);
      maxMemory *= 1024;
      maxMemory *= 1024;
      maxMemory *= 1024;

    } else if (strcmp(argv[arg], "-deleteearly") == 0) {
      deleteIntermediateEarly = true;

    } else if (strcmp(argv[arg], "-deletelate") == 0) {
      deleteIntermediateLate  = true;

    } else if (strcmp(argv[arg], "-force") == 0) {
      forceRun = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if (jobIndex == 0)
    err++;
  if (jobIdxMax == 0)
    err++;
  if (err) {
    exit(1);
  }

  //  Check if we're running or done (or crashed), then note that we're running.

  {
    char name[FILENAME_MAX];
    sprintf(name,"%s/%04d.ovs", ovlName, jobIndex);

    if ((forceRun == false) && (AS_UTL_fileExists(name, FALSE, FALSE)))
      fprintf(stderr, "Job "F_U32" is running or finished (remove '%s' or -force to try again).\n", jobIndex, name), exit(0);

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    fclose(F);	
  }

  // Get sizes of each bucket, and the final merge

  uint64   *sliceSizes    = new uint64 [fileLimit + 1];  //  For each overlap job, number of overlaps per bucket
  uint64   *bucketSizes   = new uint64 [jobIdxMax + 1];  //  For each bucket we care about, number of overlaps

  uint64    totOvl        = 0;
  uint64    numOvl        = 0;

  for (uint32 i=0; i<=jobIdxMax; i++) {
    bucketSizes[i] = 0;

    char namz[FILENAME_MAX];
    char name[FILENAME_MAX];

    sprintf(namz, "%s/bucket%04d/slice%03d.gz", ovlName, i, jobIndex);
    sprintf(name, "%s/bucket%04d/slice%03d",    ovlName, i, jobIndex);

    if ((AS_UTL_fileExists(namz, FALSE, FALSE) == false) &&
        (AS_UTL_fileExists(name, FALSE, FALSE) == false))
      //  If no file, there are no overlaps.  Skip loading the bucketSizes file.
      //  We expect the gz version to exist (that's the default in bucketizer) more frequently, so
      //  be sure to test for existence of that one first.
      continue;

    sprintf(name, "%s/bucket%04d/sliceSizes", ovlName, i);

    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "ERROR:  Failed to open %s: %s\n", name, strerror(errno)), exit(1);

    uint64 nr = AS_UTL_safeRead(F, sliceSizes, "sliceSizes", sizeof(uint64), fileLimit + 1);

    fclose(F);

    if (nr != fileLimit + 1) {
      fprintf(stderr, "ERROR: short read on '%s'.\n", name);
      fprintf(stderr, "ERROR: read "F_U64" sizes insteadof "F_U32".\n", nr, fileLimit + 1);
    }
    assert(nr == fileLimit + 1);

    fprintf(stderr, "Found "F_U64" overlaps from '%s'.\n", sliceSizes[jobIndex], name);

    bucketSizes[i] = sliceSizes[jobIndex];
    totOvl        += sliceSizes[jobIndex];
  }

  delete [] sliceSizes;
  sliceSizes = NULL;

  if (sizeof(OVSoverlap) * totOvl > maxMemory) {
    fprintf(stderr, "ERROR:  Overlaps need %.2f GB memory, but process limited (via -M) to "F_U64" GB.\n",
            sizeof(OVSoverlap) * totOvl / (1024.0 * 1024.0 * 1024.0), maxMemory >> 30);

    char name[FILENAME_MAX];
    sprintf(name,"%s/%04d.ovs", ovlName, jobIndex);

    unlink(name);

    exit(1);
  }

  fprintf(stderr, "Overlaps need %.2f GB memory, allowed to use up to (via -M) "F_U64" GB.\n",
          sizeof(OVSoverlap) * totOvl / (1024.0 * 1024.0 * 1024.0), maxMemory >> 30);

  OVSoverlap *overlapsort = new OVSoverlap [totOvl];

  //  Load all overlaps - we're guaranteed that either 'name.gz' or 'name' exists (we checked above)
  //  or funny business is happening with our files.

  for (uint32 i=0; i<=jobIdxMax; i++) {
    if (bucketSizes[i] == 0)
      continue;

    char name[FILENAME_MAX];

    sprintf(name, "%s/bucket%04d/slice%03d.gz", ovlName, i, jobIndex);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
      sprintf(name, "%s/bucket%04d/slice%03d", ovlName, i, jobIndex);

    if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
      fprintf(stderr, "ERROR: "F_U64" overlaps claim to exist in bucket '%s', but file not found.\n",
              bucketSizes[i], name);

    fprintf(stderr, "Loading "F_U64" overlaps from '%s'.\n", bucketSizes[i], name);

    BinaryOverlapFile *bof = AS_OVS_openBinaryOverlapFile(name, FALSE);
    uint64             num = 0;

    while (AS_OVS_readOverlap(bof, overlapsort + numOvl)) {
      numOvl++;
      num++;
    }

    if (num != bucketSizes[i])
      fprintf(stderr, "ERROR: expected "F_U64" overlaps, found "F_U64" overlaps.\n", bucketSizes[i], num);
    assert(num == bucketSizes[i]);

    AS_OVS_closeBinaryOverlapFile(bof);
  }

  if (numOvl != totOvl)
    fprintf(stderr, "ERROR: read "F_U64" overlaps, expected "F_U64"\n", numOvl, totOvl);
  assert(numOvl == totOvl);

  if (deleteIntermediateEarly) {
    char name[FILENAME_MAX];

    fprintf(stderr, "Removing inputs.\n");
    for (uint32 i=0; i<=jobIdxMax; i++) {
      if (bucketSizes[i] == 0)
        continue;

      sprintf(name, "%s/bucket%04d/slice%03d.gz", ovlName, i, jobIndex);
      AS_UTL_unlink(name);

      sprintf(name, "%s/bucket%04d/slice%03d", ovlName, i, jobIndex);
      AS_UTL_unlink(name);
    }
  }

  //  Sort the overlaps - at least on FreeBSD 8.2 with gcc46, the parallel STL sort
  //  algorithms are NOT inplace.  Restrict to sequential sorting.
  //
  //  This sort takes at most 2 minutes on 7gb of overlaps.
  //
  fprintf(stderr, "Sorting.\n");

#ifdef _GLIBCXX_PARALLEL
  //  If we have the parallel STL, don't use it!  Sort is not inplace!
  __gnu_sequential::sort(overlapsort, overlapsort + numOvl);
#else
  sort(overlapsort, overlapsort + numOvl);
#endif

  //  Output to store format

  fprintf(stderr, "Writing output.\n");
  writeOverlaps(ovlName, overlapsort, numOvl, jobIndex);

  delete [] overlapsort;

  if (deleteIntermediateLate) {
    char name[FILENAME_MAX];

    fprintf(stderr, "Removing inputs.\n");
    for (uint32 i=0; i<=jobIdxMax; i++) {
      if (bucketSizes[i] == 0)
        continue;

      sprintf(name, "%s/bucket%04d/slice%03d.gz", ovlName, i, jobIndex);
      AS_UTL_unlink(name);

      sprintf(name, "%s/bucket%04d/slice%03d", ovlName, i, jobIndex);
      AS_UTL_unlink(name);
    }
  }

  delete [] bucketSizes;
}
