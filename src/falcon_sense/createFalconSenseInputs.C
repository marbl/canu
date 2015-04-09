const char *mainid = "$Id:  $";

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include <vector>
#include <algorithm>

using namespace std;

int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *ovlName   = 0L;
  char             *tigName   = 0L;
  uint32            tigVers   = 0;

  uint32            errorRate = AS_OVS_encodeQuality(0.015);

  char             *outputPrefix  = NULL;

  argc = AS_configure(argc, argv);

  uint32            iidMin   = 0;
  uint32            iidMax   = UINT32_MAX;

  uint32            numReadsPer   = 0;
  uint32            numPartitions = 128;

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

  //  Open gkpStore.  Pretty much the first thing we always do.

  gkStore  *gkpStore = new gkStore(gkpName);

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

  fprintf(stderr, "Will partition "F_U32" total child reads into "F_U32" partitions.\n",
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

  //  Output falcon input.  The format is:
  //    readID fastaseq
  //    mapid1 fastaseq
  //    mapid2 fastaseq
  //    ..
  //    ++                 <- means to call consensus for 'readID' with the fastaseq supplied so far
  //    --                 <- means to stop processing, end of file

  gkRead      *read;
  gkReadData   readData;

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

    gkpStore->gkStore_loadReadData(ti, &readData);

    if (partFile[pp] == NULL) {
      char  name[FILENAME_MAX];

      sprintf(name, "%s%04d", outputPrefix, pp);  //  Sync'd with ca3g/CorrectReads.pm

      errno = 0;
      partFile[pp] = fopen(name, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
    }

    fprintf(partFile[pp], "read"F_U32" %s\n", ti, readData.gkReadData_getSequence());

    for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
      gkpStore->gkStore_loadReadData(tig->getChild(cc)->ident(), &readData);

      fprintf(partFile[pp], "data"F_U32" %s\n", tig->getChild(cc)->ident(), readData.gkReadData_getSequence());
    }

    fprintf(partFile[pp], "+ +\n");
  }


  for (uint32 pp=0; pp<numPartitions; pp++) {
    if (partFile[pp] == NULL)
      continue;

    fprintf(partFile[pp], "- -\n");
    fclose(partFile[pp]);
  }


  delete tigStore;
  delete gkpStore;

  return(0);
}
