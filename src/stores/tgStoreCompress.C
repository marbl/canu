
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"
#include "tgStore.H"


void
operationCompress(char *tigName, int tigVers) {
  tgStore    *tigStore  = new tgStore(tigName, tigVers);
  uint32      nErrors   = 0;
  uint32      nCompress = 0;

  //  Fail if this isn't the latest version.  If we try to compress something that isn't the latest
  //  version, versions after this still point to the uncompressed tigs.
  //
  //  Function never written - is this still a problem?  (18 APR 2018)


  //  Check that we aren't going to pull a tig out of the future and place it in the past.

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (tigStore->getVersion(ti) > tigVers) {
      fprintf(stderr, "WARNING:  Attempt to move future unitig " F_U32 " from version " F_U32 " to previous version %d.\n",
              ti, tigStore->getVersion(ti), tigVers);
      nErrors++;
    } else if (tigStore->getVersion(ti) < tigVers) {
      nCompress++;
    }
  }

  if (nErrors > 0) {
    fprintf(stderr, "Store can't be compressed; probably trying to compress to something that isn't the latest version.\n");
    fprintf(stderr, "  " F_U32 " tigs failed; " F_U32 " compressable\n", nErrors, nCompress);
    delete tigStore;
    exit(1);
  }


  //  Actually do the moves.

  if (nCompress > 0) {
    delete tigStore;
    tigStore = new tgStore(tigName, tigVers, tgStoreModify);
  }

  if (nCompress > 0) {
    fprintf(stderr, "Compressing " F_U32 " tigs into version %d\n", nCompress, tigVers);

    for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
      if ((ti % 1000000) == 0)
        fprintf(stderr, "tig %d\n", ti);

      if (tigStore->isDeleted(ti)) {
        continue;
      }

      if (tigStore->getVersion(ti) == tigVers)
        continue;

      tgTig *tig = tigStore->loadTig(ti);

      if (tig == NULL)
        continue;

      tigStore->insertTig(tig, true);
      tigStore->unloadTig(ti);
    }
  }

  //  Clean up the older files.

  if (nCompress > 0) {
    for (uint32 version=1; version<tigVers; version++) {
      fprintf(stderr, "Purge version " F_U32 ".\n", version);
      tigStore->purgeVersion(version);
    }
  }

  //  And the newer files.

  delete tigStore;
}






int
main (int argc, char **argv) {
  char            *seqName   = NULL;
  char            *tigName   = NULL;
  int32            tigVers   = -1;
  vector<char *>   tigInputs;
  tgStoreType      tigType   = tgStoreModify;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (seqName == NULL) || (tigName == NULL) || (tigInputs.size() == 0)) {
    fprintf(stderr, "usage: %s -S <seqStore> -T <tigStore> <v>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S <seqStore>         Path to a sequence store\n");
    fprintf(stderr, "  -T <tigStore> <v>     Path to a tigStore and version to add tigs to\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Remove store versions before <v>.  Data present in versions before <v>\n");
    fprintf(stderr, "  are copied to version <v>.  Files for the earlier versions are removed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  WARNING!  This code HAS NOT been tested with canu.\n");
    fprintf(stderr, "\n");

    if (seqName == NULL)
      fprintf(stderr, "ERROR:  no sequence store (-S) supplied.\n");
    if (tigName == NULL)
      fprintf(stderr, "ERROR:  no tig store (-T) supplied.\n");

    exit(1);
  }

  operationCompress(tigName, tigVers);

  exit(0);
}
