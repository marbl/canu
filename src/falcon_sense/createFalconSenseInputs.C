
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
 *    Brian P. Walenz from 2015-APR-09 to 2015-AUG-14
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "outputFalcon.H"

#include <vector>
#include <algorithm>

using namespace std;


int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *ovlName   = 0L;
  char             *tigName   = 0L;
  uint32            tigVers   = 0;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix  = NULL;

  argc = AS_configure(argc, argv);

  uint32            iidMin   = 0;
  uint32            iidMax   = UINT32_MAX;

  uint32            numReadsPer   = 0;
  uint32            numPartitions = 128;

  bool              trimToAlign  = true;

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax  = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-n") == 0) {
      numReadsPer   = atoi(argv[++arg]);
      numPartitions = 0;

    } else if (strcmp(argv[arg], "-p") == 0) {
      numReadsPer   = 0;
      numPartitions = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (err) {
    exit(1);
  }
  if (gkpName == NULL || tigName == NULL) {
    exit(1);
  }

  //  Open gkpStore.  Pretty much the first thing we always do.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);

  //  Open tigStore, check ranges.

  tgStore  *tigStore = new tgStore(tigName, tigVers);

  uint32   nTigs = tigStore->numTigs();

  if (nTigs <= iidMax)
    iidMax = nTigs - 1;

  //  Count how many reads are referenced in tigs from iidMin to iidMax.  These unitigs
  //  are special in that the can contain duplicate reads, so there will be many more
  //  reads referenced than the number of reads in gkpStore.

  uint64  nReadsInTigs = 0;

  for (uint32 ti=iidMin; ti<=iidMax; ti++)
    nReadsInTigs += tigStore->getNumChildren(ti);

  //  Decide how many partitions there should be.  Rather easy if the value is supplied,
  //  but if not, compute it from the number of reads per partition.

  if (numReadsPer > 0)
    numPartitions = nReadsInTigs / numReadsPer + 1;

  fprintf(stderr, "Will partition "F_U64" total child reads into "F_U32" partitions.\n",
          nReadsInTigs, numPartitions);

  //  Decide on a partitioning, based on total reads per tig.

  uint32  *tigToPart     = new uint32 [nTigs];
  uint32  *nReadsPerPart = new uint32 [numPartitions + 1];

  memset(tigToPart,     0, sizeof(uint32) * (nTigs));
  memset(nReadsPerPart, 0, sizeof(uint32) * (numPartitions + 1));

  //  Grab the number of reads per tig, again, but let us sort it to do a simple
  //  greedy partitioning.

  vector<pair<uint32,uint32> >  readsPerTig;

  for (uint32 ti=iidMin; ti<=iidMax; ti++)
    readsPerTig.push_back(pair<uint32,uint32>(tigStore->getNumChildren(ti), ti));

  sort(readsPerTig.rend(), readsPerTig.rbegin());

  //  Put the next unitig in the most empty partition.  Definitely better algorithms exist...

  for (uint32 ii=0; ii<readsPerTig.size(); ii++) {
    uint32 nReads = readsPerTig[ii].first;
    uint32 ti     = readsPerTig[ii].second;

    //  Find smallest partition

    uint32 s = 1;

    for (uint32 pp=2; pp <= numPartitions; pp++) {
      if (nReadsPerPart[pp] < nReadsPerPart[s])
        s = pp;
    }

    nReadsPerPart[s] += nReads;
    tigToPart[ti]     = s;
  }

  //  Output falcon input.

  gkRead      *read;
  gkReadData  *readData = new gkReadData;

  FILE       **partFile = new FILE * [numPartitions + 1];
  memset(partFile, 0, sizeof(FILE *) * (numPartitions + 1));

  for (uint32 ti=iidMin; ti<=iidMax; ti++) {
    tgTig *tig = tigStore->loadTig(ti);

    if (tig == NULL)
      continue;

    if (tig->numberOfChildren() == 0)
      continue;

    uint32  pp = tigToPart[ti];

    assert(pp > 0);

    if (partFile[pp] == NULL) {
      char  name[FILENAME_MAX];

      sprintf(name, "%s%04d", outputPrefix, pp);  //  Sync'd with canu/CorrectReads.pm

      errno = 0;
      partFile[pp] = fopen(name, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
    }

    outputFalcon(gkpStore, tig, trimToAlign, partFile[pp], readData);
  }

  delete readData;

  for (uint32 pp=0; pp<=numPartitions; pp++) {
    if (partFile[pp] == NULL)
      continue;

    fprintf(partFile[pp], "- -\n");
    fclose(partFile[pp]);
  }

  delete tigStore;
  delete [] partFile;
  delete [] tigToPart;
  delete [] nReadsPerPart;

  gkpStore->gkStore_close();

  return(0);
}
