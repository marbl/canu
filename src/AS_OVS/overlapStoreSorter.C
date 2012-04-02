
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

const char *mainid = "$Id: overlapStoreSorter.C,v 1.1 2012-04-02 10:58:04 brianwalenz Exp $";

#include "AS_PER_gkpStore.h"

#include "overlapStore.h"

#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"
#include "AS_OBT_acceptableOverlap.h"

#include <ctype.h>
#include <unistd.h>  //  sysconf()

#include <vector>
#include <algorithm>

using namespace std;

#define AS_OVS_CURRENT_VERSION  2

#define WITH_GZIP 1
#define DELETE_INTERMEDIATE


//  This should be private to AS_OVS
//
void
writeOverlaps(char                *ovlName,
              OVSoverlap          *overlapsort,
              uint32               numOvl,
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

  {
    sprintf(name,"%s/%04d.ovs", ovlName, jobIndex);

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);
  
    AS_UTL_safeWrite(F, &ovs, "Partition ovs file", sizeof(OverlapStoreInfo), 1);

    fclose(F);	

    fprintf(stderr, "Wrote "F_U64" overlaps into '%s'\n", ovs.numOverlapsTotal, name);
  }

  AS_OVS_closeBinaryOverlapFile(bof);    
}
















int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  char           *gkpName      = NULL;
  uint32          fileLimit    = 512;

  uint32          jobIndex     = 0;
  uint32          jobIdxMax    = 0;

  uint32          nThreads = 4;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-job") == 0) {
      jobIndex  = atoi(argv[++arg]);
      jobIdxMax = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if (gkpName == NULL)
    err++;
  if (jobIndex == 0)
    err++;
  if (jobIdxMax == 0)
    err++;
  if (err) {
    exit(1);
  }

  gkStore *gkp            = new gkStore(gkpName, FALSE, FALSE);

  // Get sizes of each bucket, and the final merge

  uint64   *dumpLength    = new uint64 [fileLimit];      //  For each overlap job, number of overlaps per bucket
  uint64   *bucketSizes   = new uint64 [jobIdxMax + 1];  //  For each bucket we care about, number of overlaps

  uint64    totOvl        = 0;
  uint64    numOvl        = 0;

  for (uint32 i=0; i<=jobIdxMax; i++) {
    bucketSizes[i] = 0;

    char name[FILENAME_MAX];
    sprintf(name, "%s/unsorted%04d/tmp.sort.%03d%s", ovlName, i, jobIndex, (WITH_GZIP) ? ".gz" : "");

    if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
      //  If no file, there are no overlaps.  Skip loading the bucketSizes file.
      continue;

    sprintf(name, "%s/unsorted%04d/bucketSizes", ovlName, i);

    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "ERROR:  Failed to open %s: %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeRead(F, dumpLength, "dumpLength", sizeof(uint64), fileLimit);

    fclose(F);

    fprintf(stderr, "Found "F_U64" overlaps from '%s'.\n", dumpLength[jobIndex], name);

    bucketSizes[i] = dumpLength[jobIndex];
    totOvl        += dumpLength[jobIndex];
  }

  delete [] dumpLength;
  dumpLength = NULL;

  OVSoverlap *overlapsort = new OVSoverlap [totOvl];

  //  Load all overlaps

  for (uint32 i=0; i<=jobIdxMax; i++) {
    if (bucketSizes[i] == 0)
      continue;

    char name[FILENAME_MAX];
    sprintf(name, "%s/unsorted%04d/tmp.sort.%03d%s", ovlName, i, jobIndex, (WITH_GZIP) ? ".gz" : "");

    fprintf(stderr, "Loading "F_U64" overlaps from '%s'.\n", bucketSizes[i], name);

    BinaryOverlapFile *bof = AS_OVS_openBinaryOverlapFile(name, FALSE);

    while (AS_OVS_readOverlap(bof, overlapsort + numOvl))
      numOvl++;

    AS_OVS_closeBinaryOverlapFile(bof);
  }

  if (numOvl != totOvl)
    fprintf(stderr, "ERROR: read "F_U64" overlaps, expected "F_U64"\n", numOvl, totOvl);
  assert(numOvl == totOvl);

  //  Sort the overlaps

  sort(overlapsort, overlapsort + numOvl);

  //  Output to store format

  writeOverlaps(ovlName, overlapsort, numOvl, jobIndex);

  delete [] overlapsort;

#ifdef DELETE_INTERMEDIATE
  for (uint32 i=0; i<=jobIdxMax; i++) {
  	if (bucketSizes[i] == 0)
  		continue;

    char name[FILENAME_MAX];
    sprintf(name, "%s/unsorted%04d/tmp.sort.%03d%s", ovlName, i, jobIndex, (WITH_GZIP) ? ".gz" : "");
  	AS_UTL_unlink(name);
  }
#endif

  delete [] bucketSizes;
}
