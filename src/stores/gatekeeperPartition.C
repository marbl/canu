
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
#include "tgStore.H"

//#include "AS_UTL_fileIO.H"

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

  //  Allocate space for the partitioning.

  uint32  *readToPart = new uint32 [numReads + 1];

  for (uint32 i=0; i<=numReads; i++)   //  All reads are in invalid
    readToPart[i] = UINT32_MAX;        //  partitions, initially.

  //  Run through all tigs and partition!

  uint32   partCount = 1;
  uint32   tigsCount = 0;
  uint32   readCount = 0;

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    //  Move to the next partition if needed

    if ((readCount + tig->numberOfChildren() >= readCountTarget) &&
        (readCount                           >  0)) {
      fprintf(stderr, "Partition %d has %d tigs and %d reads.\n",
              partCount, tigsCount, readCount);

      partCount++;
      tigsCount = 0;
      readCount = 0;
    }

    //  Assign all the reads in this tig to this partition.

    readCount += tig->numberOfChildren();

    for (uint32 ci=0; ci<tig->numberOfChildren(); ci++)
      readToPart[tig->getChild(ci)->ident()] = partCount;

    tigStore->unloadTig(ti);
  }

  if (readCount > 0)
    fprintf(stderr, "Partition %d has %d tigs and %d reads.\n",
            partCount, tigsCount, readCount);

  delete tigStore;

  return(readToPart);
}



int
main(int argc, char **argv) {
  char   *gkpStorePath      = NULL;
  char   *tigStorePath      = NULL;
  uint32  tigStoreVers      = 0;
  uint32  readCountTarget   = 2500;   //  No partition smaller than this
  uint32  partCountTarget   = 200;    //  No more than this many partitions

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigStorePath = argv[++arg];
      tigStoreVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      readCountTarget = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      partCountTarget = atoi(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR: unknown option '%s'\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (gkpStorePath == NULL)  err.push_back("ERROR: no gkpStore (-G) supplied.\n");
  if (tigStorePath == NULL)  err.push_back("ERROR: no partition input (-P) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G <gkpStore> -T <tigStore> <v>\n", argv[0]);
    fprintf(stderr, "  -G <gkpStore>       path to gatekeeper store\n");
    fprintf(stderr, "  -T <tigStore> <v>   path to tig store and version to be partitioned\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b <nReads>         minimum number of reads per partition (50000)\n");
    fprintf(stderr, "  -p <nPartitions>    number of partitions (200)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Create a partitioned copy of <gkpStore> and place it in <tigStore>/partitionedReads.gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "NOTE:  Path handling in this is probably quite brittle.  Due to an implementation\n");
    fprintf(stderr, "       detail, the new store must have symlinks back to the original store.  Canu \n");
    fprintf(stderr, "       wants to use relative paths, and this program tries to adjust <gkpStore> to be\n");
    fprintf(stderr, "       relative to <tigStore/partitionedReads.gkpStore.  If it fails to do this correctly,\n");
    fprintf(stderr, "       one of two (seen so far) errors will occur:\n");
    fprintf(stderr, "         Original file '.../partitionedReads.gkpStore/info' doesn't exist, won't make a link to nothing.\n");
    fprintf(stderr, "         Couldn't open '.../partitionedReads.gkpStore/libraries' for mmap: No such file or directoryn");
    fprintf(stderr, "       In both cases, try to simplify <tigStore> -- in particular, remove any '..' or '.' components -- or\n");
    fprintf(stderr, "       run this from a higher/lower directory.\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  char    gkpSourcePath[FILENAME_MAX] = {0};
  char    gkpClonePath[FILENAME_MAX]  = {0};

  //  We're making a clone of the master gkpStore in the tigStore directory.

  snprintf(gkpClonePath, FILENAME_MAX, "%s/partitionedReads.gkpStore", tigStorePath);

  //  The path to the gkpStore that we want to use in the link is a wee-bit more complicated.
  //  If it's an absolute path, there's nothing we need to do.

  if (gkpStorePath[0] == '/') {
    strcpy(gkpSourcePath, gkpStorePath);
  }

  //  But if it's a relative path, we need to add a bunch of dots.  One pair to account
  //  for the directory we added above, and then more dots for each component in tigStorePath.

  else {
    char    t[FILENAME_MAX];                  //  Copy command line tigStorePath to a
    char   *p = t;                            //  local, and modifiable, space.

    strcpy(p, tigStorePath);

    strcat(gkpSourcePath, "../");             //  One for the directory we created above

    while ((p[0] != '.') || (p[1] != 0)) {    //  Many for each component in the tigStorePath.
      strcat(gkpSourcePath, "../");
      p = dirname(p);
    }

    if ((gkpStorePath[0] == '.') && (gkpStorePath[1] == '/'))   //  Finally, append the supplied
      strcat(gkpSourcePath, gkpStorePath + 2);                  //  gkpStorePath, possibly
    else                                                        //  stripping off any ./ at the
      strcat(gkpSourcePath, gkpStorePath);                      //  start.
  }

  //  Make the clone.

  gkStore::gkStore_clone(gkpSourcePath, gkpClonePath);

  //  Open the clone.

  gkStore *gkpStore = gkStore::gkStore_open(gkpClonePath, gkStore_readOnly);

  //  Scan all the tigs to build a map from read to partition.

  uint32   *partition = buildPartition(tigStorePath, tigStoreVers,
                                       readCountTarget,
                                       partCountTarget,
                                       gkpStore->gkStore_getNumReads());

  //  Dump the partition data to the store, let it build partitions.

  gkpStore->gkStore_buildPartitions(partition);

  //  That's all folks.

  delete [] partition;

  gkpStore->gkStore_close();

  exit(0);
}
