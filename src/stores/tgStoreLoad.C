
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
 *    Brian P. Walenz from 2015-AUG-07 to 2015-AUG-14
 *      are Copyright 2015 Battelle National Biodefense Institute, and
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



void
operationBuild(char   *buildName,
               char   *tigName,
               uint32  tigVers) {

  errno = 0;
  FILE *F = fopen(buildName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", buildName, strerror(errno)), exit(1);

  if (AS_UTL_fileExists(tigName, TRUE, TRUE)) {
    fprintf(stderr, "ERROR: '%s' exists, and I will not clobber an existing store.\n", tigName);
    exit(1);
  }

  tgStore *tigStore = new tgStore(tigName);
  tgTig    *tig      = new tgTig();

  for (int32 v=1; v<tigVers; v++)
    tigStore->nextVersion();

  while (tig->loadLayout(F) == true) {
    if (tig->numberOfChildren() == 0)
      continue;

    //  The log isn't correct.  For new tigs (all of these are) we don't know the
    //  id until after it is added.  Further, if these come with id's already set,
    //  they can't be added to a new store -- they don't exist.

#if 0
    fprintf(stderr, "INSERTING tig %d (%d children) (originally ID %d)\n",
            tig->tigID(), tig->numberOfChildren(), oID);
#endif

    tigStore->insertTig(tig, false);
  }

  fclose(F);

  delete tig;
  delete tigStore;
}





int
main (int argc, char **argv) {
  char            *gkpName   = NULL;
  char            *tigName   = NULL;
  int32            tigVers   = -1;
  vector<char *>   tigInputs;
  tgStoreType      tigType   = tgStoreModify;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      AS_UTL_loadFileList(argv[++arg], tigInputs);

    } else if (strcmp(argv[arg], "-n") == 0) {
      tigType = tgStoreReadOnly;

    } else if (AS_UTL_fileExists(argv[arg])) {
      tigInputs.push_back(argv[arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL) || (tigInputs.size() == 0)) {
    fprintf(stderr, "usage: %s -G <gkpStore> -T <tigStore> <v> [input.cns]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <gkpStore>         Path to the gatekeeper store\n");
    fprintf(stderr, "  -T <tigStore> <v>     Path to the tigStore and version to add tigs to\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L <file-of-files>    Load the tig(s) from files listed in 'file-of-files'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -n                    Don't replace, just report what would have happened\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  The primary operation is to replace tigs in the store with ones in a set of input files.\n");
    fprintf(stderr, "  The input files can be either supplied directly on the command line or listed in\n");
    fprintf(stderr, "  a text file (-L).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  A new store is created if one doesn't exist, otherwise, whatever tigs are there are\n");
    fprintf(stderr, "  replaced with those in the -R file.  If version 'v' doesn't exist, it is created.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Even if -n is supplied, a new store is created if one doesn't exist.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  To add a new tig, give it a tig id of -1.  New tigs must be added to the latest version.\n");
    fprintf(stderr, "  To delete a tig, remove all children, and set the number of them to zero.\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR:  no gatekeeper store (-G) supplied.\n");
    if (tigName == NULL)
      fprintf(stderr, "ERROR:  no tig store (-T) supplied.\n");
    if (tigInputs.size() == 0)
      fprintf(stderr, "ERROR:  no input tigs (-R) supplied.\n");

    exit(1);
  }

  //  If the store doesn't exist, create one, and make a bunch of versions
  if (AS_UTL_fileExists(tigName, true, false) == false) {
    fprintf(stderr, "Creating tig store '%s' version %d\n", tigName, tigVers);

    tgStore *tigStore = new tgStore(tigName);

    for (int32 vv=1; vv<tigVers; vv++)
      tigStore->nextVersion();

    delete tigStore;
  }

  gkStore *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore *tigStore = new tgStore(tigName, tigVers, tigType);
  tgTig   *tig      = new tgTig;

  for (uint32 ff=0; ff<tigInputs.size(); ff++) {
    errno = 0;
    FILE *TI = fopen(tigInputs[ff], "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", tigInputs[ff], strerror(errno)), exit(1);

    fprintf(stderr, "Reading layouts from '%s'.\n", tigInputs[ff]);

    while (tig->loadFromStreamOrLayout(TI) == true) {

      //  Handle insertion.

      if (tig->numberOfChildren() > 0) {
        //fprintf(stderr, "INSERTING tig %d\n", tig->tigID());
        tigStore->insertTig(tig, false);
        continue;
      }

      //  Deleted already?

      if (tigStore->isDeleted(tig->tigID()) == true) {
        //fprintf(stderr, "DELETING tig %d -- ALREADY DELETED\n", tig->tigID());
        continue;
      }

      //  Really delete it then.

      //fprintf(stderr, "DELETING tig %d\n", tig->tigID());
      tigStore->deleteTig(tig->tigID());
    }

    fclose(TI);

    fprintf(stderr, "Reading layouts from '%s' completed.\n", tigInputs[ff]);
  }

  delete tig;
  delete tigStore;

  gkpStore->gkStore_close();

  exit(0);
}
