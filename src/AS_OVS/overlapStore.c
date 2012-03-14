
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

const char *mainid = "$Id: overlapStore.c,v 1.34 2012-03-14 14:58:35 gesims Exp $";

#include "overlapStore.h"
#include "AS_OVS_overlap.h"   //  Just to know the sizes of structs
#include "AS_PER_gkpStore.h"  //  Just to know clear region labels

#include <ctype.h>
#include <unistd.h>  //  sysconf()

#include <vector>

using namespace std;

int
main(int argc, char **argv) {
  uint32          operation   = OP_NONE;
  char           *storeName   = NULL;
  char           *gkpName     = NULL;
  uint32          clearRegion = AS_READ_CLEAR_ERROR;

  uint32          dumpBinary  = FALSE;
  uint32 optimizedStoreCreate = FALSE;
  uint32 optimizedOvlBucketize= FALSE;
  uint32 optimizedMergeBuckets= FALSE;
  uint32 optimizedSortMergedBuckets = FALSE;
  uint32 optimizedBuildIndexFile = FALSE;
  uint32 workingOvlDumpIndex  = 0;
  uint32 workingBucketIndex  = 0;
  double          dumpERate   = 100.0;
  float		  quality = 0.04;
  uint32          dumpType    = 0;

  uint32          bgnIID      = 0;
  uint32          endIID      = UINT32_MAX;
  uint32          qryIID      = 0;

  uint32          gs_ovlLimit = 100;
  uint32          gs_winSize  = 100;
  uint32          gs_minOvl   = 40;
  uint32          gs_into     = 5;

  uint64          memoryLimit = 512 * 1024 * 1024;
  uint32          fileLimit   = 0;
  uint32          nThreads    = 4;
  uint32          doFilterOBT = 0;
  vector<char *>  fileList;
  Ovl_Skip_Type_t ovlSkipOpt  = ALL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-c") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_BUILD;

    } else if (strcmp(argv[arg], "-m") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_MERGE;

    } else if (strcmp(argv[arg], "-d") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_DUMP;

    } else if (strcmp(argv[arg], "-p") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      qryIID      = atoi(argv[++arg]);
      storeName   = argv[++arg];
      gkpName     = argv[++arg];
      clearRegion = gkStore_decodeClearRegionLabel(argv[++arg]);
      operation   = OP_DUMP_PICTURE;

    } else if (strcmp(argv[arg], "-G") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName    = argv[++arg];
      gkpName      = argv[++arg];
      gs_ovlLimit  = atoi(argv[++arg]);
      operation    = OP_GENOME_LENGTH;

    } else if (strcmp(argv[arg], "-u") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_UPDATE_ERATES;

    } else if (strcmp(argv[arg], "-t") == 0) {
      nThreads    = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      dumpERate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-d5") == 0) {
      dumpType |= DUMP_5p;

    } else if (strcmp(argv[arg], "-d3") == 0) {
      dumpType |= DUMP_3p;

    } else if (strcmp(argv[arg], "-dC") == 0) {
      dumpType |= DUMP_CONTAINS;

    } else if (strcmp(argv[arg], "-dc") == 0) {
      dumpType |= DUMP_CONTAINED;

    } else if (strcmp(argv[arg], "-B") == 0) {
      dumpBinary = TRUE;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-q") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      bgnIID    = atoi(argv[++arg]);
      endIID    = bgnIID;
      qryIID    = atoi(argv[++arg]);
      storeName = argv[++arg];
      operation = OP_DUMP;

    } else if (strcmp(argv[arg], "-O") == 0) {
      doFilterOBT++;

    } else if (strcmp(argv[arg], "-M") == 0) {
      memoryLimit  = atoi(argv[++arg]);  //  convert first, then multiply so we don't
      memoryLimit *= 1024 * 1024;        //  overflow whatever type atoi() is.
      fileLimit    = 0;

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit    = atoi(argv[++arg]);
      memoryLimit  = 0;

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-U") == 0) {
      optimizedOvlBucketize = TRUE;
      workingOvlDumpIndex = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-Q") == 0) {
      quality = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-W") == 0) {
      optimizedSortMergedBuckets = TRUE;
      workingBucketIndex = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-P") == 0) {
      optimizedStoreCreate = TRUE;
      
    } else if (strcmp(argv[arg], "-I") == 0) {
      optimizedBuildIndexFile = TRUE;

    } else if (strcmp(argv[arg], "-L") == 0) {
      char *line;

      //  The next arg is a file with the list of files to use
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno)), exit(1);

      line = (char *)safe_malloc(sizeof(char) * FILENAME_MAX);
      fgets(line, FILENAME_MAX, F);
      while (!feof(F)) {
        chomp(line);
        fileList.push_back(line);
       line = (char *)safe_malloc(sizeof(char) * FILENAME_MAX);
        fgets(line, FILENAME_MAX, F);
      }
      safe_free(line);
      fclose(F);
    
    } else if (strcmp(argv[arg], "-i") == 0) {      
      switch (atoi(argv[++arg])) {
         case 0:
            ovlSkipOpt = NONE;
            break;
         case 1:
            ovlSkipOpt = ALL;
            break;
         case 2:
            ovlSkipOpt = INTERNAL;
            break;
         default:
            fprintf(stderr, "%s: unknown overlap ignore option '%s'. Must be one of 0, 1, or 2.\n", argv[0], argv[arg]);
            err++;
            break;
      }
    
    } else if ((argv[arg][0] == '-') && (argv[arg][1] != 0)) {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;

    } else {
      //  Assume it's an input file
      fileList.push_back(argv[arg]);
    }
    arg++;
  }
  if ((operation == OP_NONE) || (storeName == NULL) || (err)) {
    fprintf(stderr, "usage: %s -c storeName [-M x (MB) | -F files] [-t threads] [-g gkpStore] [-L list-of-ovl-files] ovl-file ...\n", argv[0]);
    fprintf(stderr, "       %s -m storeName mergeName\n", argv[0]);
    fprintf(stderr, "       %s -d storeName [-B] [-E erate] [-b beginIID] [-e endIID]\n", argv[0]);
    fprintf(stderr, "       %s -q aiid biid storeName\n", argv[0]);
    fprintf(stderr, "       %s -p iid storeName gkpStore clr\n", argv[0]);
    fprintf(stderr, "       %s -G storeName gkpStore ovlLimit\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "There are six modes of operation, selected by the first option:\n");
    fprintf(stderr, "  -c  create a new store, fails if the store exists\n");
    fprintf(stderr, "  -m  merge store mergeName into store storeName\n");
    fprintf(stderr, "  -d  dump a store\n");
    fprintf(stderr, "  -q  report the a,b overlap, if it exists.\n");
    fprintf(stderr, "  -p  dump a picture of overlaps to fragment 'iid', using clear region 'clr'.\n");
    fprintf(stderr, "  -G  estimate genome length\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CREATION - create a new store from raw overlap files\n");
    fprintf(stderr, "  -P           (alpha) Optimized overlapstore create.\n");
    fprintf(stderr, "  -U x         (alpha) Optimized grid enabled bucketizer  (bucketize xth overlap file)\n");
    fprintf(stderr, "          -Q f         Specify error in overlap filtering for above: Default 0.04\n");
    fprintf(stderr, "  -W x         (alpha) Sort buckets created with -U and dump to storefile.\n");
    fprintf(stderr, "  -I           (alpha) Create index of all storefile buckets.\n");
    fprintf(stderr, "  -O           Filter overlaps for OBT.\n");
    fprintf(stderr, "  -M x         Use 'x'MB memory for sorting overlaps (conflicts with -F).\n");
    fprintf(stderr, "  -F x         Use 'x' files for sorting overlaps (conflicts with -M).\n");
    fprintf(stderr, "  -t t         Use 't' threads for sorting overlaps.\n");
    fprintf(stderr, "  -L f         Read overlaps from files listed in 'f'.\n");
    fprintf(stderr, "  -i x         Ignore overlaps to closure reads; x is:\n");
    fprintf(stderr, "                 0 Delete no overlaps.\n");
    fprintf(stderr, "                 1 Delete all overlaps to closure read (default).\n");
    fprintf(stderr, "                 2 Delete only overlaps between closure reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "MERGING - merge two stores into one\n");
    fprintf(stderr, "  -m storeName mergeName   Merge the store 'mergeName' into 'storeName'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING - report overlaps in the store\n");
    fprintf(stderr, "  -B                Dump the store as binary, suitable for input to create a new store.\n");
    fprintf(stderr, "  -E erate          Dump only overlaps <= erate error.\n");
    fprintf(stderr, "  -d5               Dump only overlaps off the 5' end of the A frag.\n");
    fprintf(stderr, "  -d3               Dump only overlaps off the 3' end of the A frag.\n");
    fprintf(stderr, "  -dC               Dump only overlaps that are contained in the A frag (B contained in A).\n");
    fprintf(stderr, "  -dc               Dump only overlaps that are containing the A frag (A contained in B).\n");
    fprintf(stderr, "  -b beginIID       Start dumping at 'beginIID'.\n");
    fprintf(stderr, "  -e endIID         Stop dumping after 'endIID'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "QUERYING - quickly ask if an overlap exists\n");
    fprintf(stderr, "  -q aiid biid storeName\n");
    fprintf(stderr, "                    If an overlap between fragments 'aiid' and 'biid' exists, it is printed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING PICTURES - draw a multi-alignment-like picture for a single fragment and its overlaps\n");
    fprintf(stderr, "  -p iid storeName gkpStore clr\n");
    fprintf(stderr, "                    clr is usually OBTINITIAL for obtStore.\n");
    fprintf(stderr, "                    clr is usually OBTCHIMERA for ovlStore when OBT is used.\n");
    fprintf(stderr, "                    clr is usually CLR        for ovlStore when OBT is not used.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Limited to %d open files.\n", (int)sysconf(_SC_OPEN_MAX));
    fprintf(stderr, "\n");
    fprintf(stderr, "OVSoverlap     %d bytes\n", (int)sizeof(OVSoverlap));
    fprintf(stderr, "OVSoverlapINT  %d bytes\n", (int)sizeof(OVSoverlapINT));
    fprintf(stderr, "OVSoverlapDAT  %d bytes\n", (int)sizeof(OVSoverlapDAT));
    fprintf(stderr, "               %d %d %d\n", (int)sizeof(struct OVSoverlapOVL), (int)sizeof(struct OVSoverlapMER), (int)sizeof(struct OVSoverlapOBT));
    fprintf(stderr, "AS_OVS_NWORDS  %d\n", AS_OVS_NWORDS);
    exit(1);
  }
  if ((fileList.size() == 0) && (operation == OP_BUILD)) {
    fprintf(stderr, "No input files?\n");
    exit(1);
  }
  if (dumpType == 0)
    dumpType = DUMP_5p | DUMP_3p | DUMP_CONTAINED | DUMP_CONTAINS;

  switch (operation) {
    case OP_BUILD:
      if (optimizedOvlBucketize) {
	   fprintf(stderr,"Using optimized grid enabled bucketizer\n");
           BucketizeOvlGES(storeName, gkpName, memoryLimit, fileLimit, nThreads, doFilterOBT, fileList.size(), &fileList[0],quality,  ovlSkipOpt, workingOvlDumpIndex );
      } else if (optimizedSortMergedBuckets) {
          sortDistributedBucketGES(storeName, gkpName, memoryLimit, fileLimit,nThreads, fileList.size(), &fileList[0], workingBucketIndex );
      } else if (optimizedBuildIndexFile) {
          buildStoreIndexGES2(storeName, gkpName, memoryLimit, fileLimit, nThreads, fileList.size(), &fileList[0]);	  
      } else if (optimizedStoreCreate) { 
	   fprintf(stderr,"Using optimized store build\n");
           buildStoreGES(storeName, gkpName, memoryLimit, fileLimit, nThreads, doFilterOBT, fileList.size(), &fileList[0], ovlSkipOpt);
      } else {	
           buildStore(storeName, gkpName, memoryLimit, fileLimit, nThreads, doFilterOBT, fileList.size(), &fileList[0], ovlSkipOpt);
      }
      break;
    case OP_MERGE:
      mergeStore(storeName, fileList[0]);
      break;
    case OP_DUMP:
      dumpStore(storeName, dumpBinary, dumpERate, dumpType, bgnIID, endIID, qryIID);
      break;
    case OP_DUMP_PICTURE:
      dumpPicture(storeName, gkpName, clearRegion, dumpERate, dumpType, qryIID);
      break;
    case OP_GENOME_LENGTH:
      estimateGenomeLength(storeName,
                           gkpName,
                           gs_ovlLimit,
                           bgnIID,
                           endIID,
                           gs_into,
                           gs_winSize,
                           gs_minOvl);
      break;
    case OP_UPDATE_ERATES:
      updateErates(storeName, fileList[0]);
      break;
    default:
      break;
  }

  exit(0);
}
