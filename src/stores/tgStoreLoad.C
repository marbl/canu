
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
dumpFile(vector<char const *>  &tigInputs) {
  tgTig   *tig      = new tgTig;

  for (uint32 ff=0; ff<tigInputs.size(); ff++) {
    fprintf(stderr, "Reading layouts from '%s'.\n", tigInputs[ff]);

    FILE *TI = AS_UTL_openInputFile(tigInputs[ff]);

    while (tig->loadFromStreamOrLayout(TI) == true)
      tig->dumpLayout(stdout);

    AS_UTL_closeFile(TI, tigInputs[ff]);

    fprintf(stderr, "Reading layouts from '%s' completed.\n", tigInputs[ff]);
  }

  delete tig;
}



void
testFile(vector<char const *>  &tigInputs) {
  tgTig   *tig      = new tgTig;

  for (uint32 ff=0; ff<tigInputs.size(); ff++) {
    FILE *TI = AS_UTL_openInputFile(tigInputs[ff]);

    while (tig->loadFromStreamOrLayout(TI) == true) {
      int32  minP = INT32_MAX;
      int32  maxP = 0;

      for (uint32 ii=0; ii<tig->numberOfChildren(); ii++) {
        minP = min(minP, tig->getChild(ii)->min());
        maxP = min(maxP, tig->getChild(ii)->max());
      }

      if ((minP != 0) ||
          (maxP != tig->length()))
        fprintf(stdout, "BAD  %8u  layout %8d %8d length %8u\n",
                tig->tigID(), minP, maxP, tig->length());
      else
        fprintf(stdout, "GOOD %8u\n", tig->tigID());
    }

    AS_UTL_closeFile(TI, tigInputs[ff]);
  }

  delete tig;
}



void                                 //  Create a new store if one doesn't exist, then
createStore(char const *tigName,     //  add empty versions until we get to the one
            int32       tigVers) {   //  we're supposed to load into.

  if (directoryExists(tigName) == false) {
    fprintf(stderr, "Creating tig store '%s' version %d\n", tigName, tigVers);

    tgStore *tigStore = new tgStore(tigName);

    for (int32 vv=1; vv<tigVers; vv++)
      tigStore->nextVersion();

    delete tigStore;
  }
}



void
loadTigs(char const            *seqName,
         char const            *tigName,
         int32                  tigVers,
         tgStoreType            tigType,
         vector<char const *>  &tigInputs) {
  sqStore *seqStore = new sqStore(seqName);
  tgStore *tigStore = new tgStore(tigName, tigVers, tigType);
  tgTig   *tig      = new tgTig;

  for (uint32 ff=0; ff<tigInputs.size(); ff++) {
    fprintf(stderr, "Reading layouts from '%s'.\n", tigInputs[ff]);

    FILE *TI = AS_UTL_openInputFile(tigInputs[ff]);

    while (tig->loadFromStreamOrLayout(TI) == true) {
      if (tig->numberOfChildren() > 0)                       //  Insert it!
        tigStore->insertTig(tig, false);

      else if (tigStore->isDeleted(tig->tigID()) == false)   //  Delete it!
        tigStore->deleteTig(tig->tigID());
    }

    AS_UTL_closeFile(TI, tigInputs[ff]);

    fprintf(stderr, "Reading layouts from '%s' completed.\n", tigInputs[ff]);
  }

  delete tig;
  delete tigStore;

  delete seqStore;
}



int
main (int argc, char **argv) {
  char const           *seqName       = NULL;
  char const           *tigName       = NULL;
  int32                 tigVers       = -1;
  vector<char const *>  tigInputs;
  char const           *tigInputsFile = NULL;
  tgStoreType           tigType       = tgStoreModify;
  bool                  doDumpFile    = false;
  bool                  doTestFile    = false;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      tigInputsFile = argv[++arg];
      AS_UTL_loadFileList(tigInputsFile, tigInputs);

    } else if (strcmp(argv[arg], "-n") == 0) {
      tigType = tgStoreReadOnly;

    } else if (strcmp(argv[arg], "-dump") == 0) {
      doDumpFile = true;

    } else if (strcmp(argv[arg], "-test") == 0) {
      doTestFile = true;

    } else if (fileExists(argv[arg])) {
      tigInputs.push_back(argv[arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((doDumpFile == false) && (doTestFile == false) && (seqName == NULL))
    err.push_back("ERROR:  no sequence store (-S) supplied.\n");
  if ((doDumpFile == false) && (doTestFile == false) && (tigName == NULL))
    err.push_back("ERROR:  no tig store (-T) supplied.\n");
  if (tigInputs.size() == 0)
    err.push_back("ERROR:  no input tig files supplied on command line or via -L option.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S <seqStore> -T <tigStore> <v> [input.cns]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S <seqStore>         Path to the sequence store\n");
    fprintf(stderr, "  -T <tigStore> <v>     Path to the tigStore and version to add tigs to\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L <file-of-files>    Load the tig(s) from files listed in 'file-of-files'\n");
    fprintf(stderr, "                        (WARNING: program will succeed if this file is empty)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -n                    Don't load into store, just report what would have happened\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -dump                 Dump the cns files as ASCII, don't load into store\n");
    fprintf(stderr, "  -test                 Test the cns files for various errors, don't load into store\n");
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

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }


  if (doDumpFile == true) {
    dumpFile(tigInputs);
  }

  else if (doTestFile == true) {
    testFile(tigInputs);
  }

  else {
    createStore(tigName, tigVers);
    loadTigs(seqName, tigName, tigVers, tigType, tigInputs);
  }

  exit(0);
}
