
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
 *    Brian P. Walenz from 2014-DEC-23 to 2015-MAR-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "AS_UTL_fileIO.H"


int
main(int argc, char **argv) {
  char   *gkpStoreName      = NULL;
  char   *partitionFile     = NULL;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-P") == 0) {
      partitionFile = argv[++arg];

    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (partitionFile == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -P partitionMapFile\n", argv[0]);
    fprintf(stderr, "  -G gkpStore         path to gatekeeper store\n");
    fprintf(stderr, "  -P partFile         file mapping read ID to partiton\n");
    fprintf(stderr, "                      format: 'partition readID'\n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-G) supplied.\n");
    if (partitionFile == NULL)
      fprintf(stderr, "ERROR: no partition input (-P) supplied.\n");
    exit(1);
  }

  //  Open a READ ONLY store.  This prevents us from mucking with the non-partitioned reads
  //  (like, by changing the read ID or pointer to the blob).  We don't need it opened
  //  for writing anyway.

  gkStore    *gkpStore  = gkStore::gkStore_open(gkpStoreName, gkStore_readOnly);
  uint32      numReads  = gkpStore->gkStore_getNumReads();

  uint32     *partition = new uint32 [numReads + 1];

  //  Set all partitions to invalid.

  for (uint32 i=0; i<=numReads; i++)
    partition[i] = UINT32_MAX;

  //  Read the partition file

  errno = 0;
  FILE *F = fopen(partitionFile, "r");
  if (errno)
    fprintf(stderr, "GKP Error: Build_Partition()-- failed to open '%s': %s\n", partitionFile, strerror(errno)), exit(1);

  while (!feof(F)) {
    uint32  i, p;

    if (2 == fscanf(F, " "F_U32" "F_U32" ", &p, &i)) {
      assert(i            <= numReads);
      assert(partition[i] == UINT32_MAX);

      partition[i] = p;
    }
  }
  fclose(F);

  //  Dump the partition data to the store, let it build partitions.

  gkpStore->gkStore_buildPartitions(partition);

  //  That's all folks.

  delete [] partition;

  gkpStore->gkStore_close();

  exit(0);
}
