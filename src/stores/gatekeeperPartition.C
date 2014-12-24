
const char *mainid = "$Id:  $";

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
    fprintf(stderr, "  -P gkpStore         file mapping read ID to partiton\n");
    fprintf(stderr, "                      format: 'readID partition'\n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-G) supplied.\n");
    if (partitionFile == NULL)
      fprintf(stderr, "ERROR: no partition input (-P) supplied.\n");
    exit(1);
  }

  gkStore    *gkpStore  = new gkStore(gkpStoreName, gkStore_modify);  //  Open for modify!
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

    if (2 == fscanf(F, " "F_U32" "F_U32" ", &i, &p))
      partition[i] = p;
  }
  fclose(F);

  //  Dump the partition data to the store, let it build partitions.

  gkpStore->gkStore_buildPartitions(partition);

  //  That's all folks.

  delete [] partition;
  delete    gkpStore;

  exit(0);
}
