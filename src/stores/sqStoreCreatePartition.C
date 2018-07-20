
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
 *    src/stores/gatekeeperPartition.C
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

#include "sqStore.H"
#include "tgStore.H"

//#include "files.H"

#include <libgen.h>


uint32 *
buildPartition(char    *tigStoreName,
               uint32   tigStoreVers,
               uint32   readCountTarget,
               uint32   partCountTarget,
               uint32   numReads) {
  tgStore *tigStore   = new tgStore(tigStoreName, tigStoreVers);

  //  Decide on how many reads per partition.  We take two targets, the partCountTarget
  //  is used to decide how many partitions to make, but if there are too few reads in
  //  each partition, we'll reset to readCountTarget.

  if (readCountTarget < numReads / partCountTarget)
    readCountTarget = numReads / partCountTarget;

  //  Figure out how many partitions we'll make, then spread the reads equally through them.

  uint32  numParts = (uint32)ceil((double)numReads / readCountTarget);

  readCountTarget = 1 + numReads / numParts;

  fprintf(stderr, "For %u reads, will make %u partition%s with up to %u reads%s.\n",
          numReads,
          (numParts),
          (numParts == 1) ? "" : "s",
          readCountTarget,
          (numParts == 1) ? "" : " in each");
  fprintf(stderr, "\n");

  //  Allocate space for the partitioning.

  uint32  *readToPart = new uint32 [numReads + 1];

  for (uint32 i=0; i<=numReads; i++)   //  All reads are in invalid
    readToPart[i] = UINT32_MAX;        //  partitions, initially.

  //  Run through all tigs and partition!

  uint32   partCount  = 1;
  uint32   tigsCount  = 0;
  uint32   readCount  = 0;
  uint32   longest    = 0;

  uint32   totalTigs  = 0;
  uint32   totalReads = 0;
  uint32   longestG   = 0;   //  Globally longest

  fprintf(stderr, "Partition      Tigs     Reads   Longest\n");
  fprintf(stderr, "--------- --------- --------- ---------\n");

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    //  Move to the next partition if needed

    if ((readCount + tig->numberOfChildren() >= readCountTarget) &&
        (readCount                           >  0)) {
      fprintf(stderr, "%9u %9u %9u %9u\n", partCount, tigsCount, readCount, longest);

      partCount++;
      tigsCount = 0;
      readCount = 0;
      longest   = 0;
    }

    //  Assign all the reads in this tig to this partition.

    tigsCount  += 1;
    readCount  += tig->numberOfChildren();

    totalTigs  += 1;
    totalReads += tig->numberOfChildren();

    longest  = max(longest,  tig->length());
    longestG = max(longestG, tig->length());

    //if (longest < tig->length())
    //  longest = tig->length();

    for (uint32 ci=0; ci<tig->numberOfChildren(); ci++)
      readToPart[tig->getChild(ci)->ident()] = partCount;

    tigStore->unloadTig(ti);
  }

  if (readCount > 0)
    fprintf(stderr, "%9u %9u %9u %9u\n", partCount, tigsCount, readCount, longest);

  fprintf(stderr, "--------- --------- --------- ---------\n");
  fprintf(stderr, "          %9u %9u %9u (partitioned)\n", totalTigs, totalReads, longestG);
  fprintf(stderr, "                    %9u           (unpartitioned)\n", numReads - totalReads);
  fprintf(stderr, "\n");

  delete tigStore;

  return(readToPart);
}



int
main(int argc, char **argv) {
  char     *seqStorePath                = NULL;
  char      seqClonePath[FILENAME_MAX]  = { 0 };
  char     *tigStorePath                = NULL;
  uint32    tigStoreVers                = 0;
  uint32    readCountTarget             = 2500;   //  No partition smaller than this
  uint32    partCountTarget             = 200;    //  No more than this many partitions
  bool      doDelete                    = false;

  sqStore  *seqStore                    = NULL;
  uint32   *partition                   = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigStorePath = argv[++arg];
      tigStoreVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      readCountTarget = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      partCountTarget = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      tigStorePath = argv[++arg];
      tigStoreVers = 1;
      doDelete = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR: unknown option '%s'\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((seqStorePath == NULL) &&
      (doDelete == false))       err.push_back("ERROR: no seqStore (-S) supplied.\n");
  if (tigStorePath == NULL)      err.push_back("ERROR: no tigStore (-T) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [-S <seqStore> -T <tigStore> <v>] ...\n", argv[0]);
    fprintf(stderr, "       %s [-D <tigStore>]\n", argv[0]);
    fprintf(stderr, "  -S <seqStore>       path to sequence store\n");
    fprintf(stderr, "  -T <tigStore> <v>   path to tig store and version to be partitioned\n");
    fprintf(stderr, "  -D <tigStore>       remove a partitioned seqStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b <nReads>         minimum number of reads per partition (50000)\n");
    fprintf(stderr, "  -p <nPartitions>    number of partitions (200)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Create a partitioned copy of <seqStore> and place it in <tigStore>/partitionedReads.seqStore\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }


  snprintf(seqClonePath, FILENAME_MAX, "%s/partitionedReads.seqStore", tigStorePath);

  //  If deleting, delete.

  if (doDelete == true) {
    seqStore = sqStore::sqStore_open(seqClonePath, sqStore_readOnly);
    seqStore->sqStore_deletePartitions();
  }

  //  Otherwise, partitioning, so partition.

  else {
    seqStore = sqStore::sqStore_open(seqStorePath,                       //  Open the store, preparing it for
                                     seqClonePath);                      //  a copy to the partitioned version.

    partition = buildPartition(tigStorePath, tigStoreVers,               //  Scan all the tigs
                               readCountTarget,                          //  to build a map from
                               partCountTarget,                          //  read to partition.
                               seqStore->sqStore_getNumReads());

    seqStore->sqStore_buildPartitions(partition);                        //  Build partitions.
  }

  //  Cleanp and bye.

  delete [] partition;

  seqStore->sqStore_close();

  fprintf(stderr, "Bye.\n");
  exit(0);
}
