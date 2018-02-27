
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
 *    Brian P. Walenz beginning on 2017-OCT-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"



int
main (int argc, char **argv) {
  char            *gkpName       = NULL;
  char            *corName       = NULL;
  int32            corVers       = 1;

  vector<char *>   corInputs;
  char            *corInputsFile = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corName = argv[++arg];
      //corVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      corInputsFile = argv[++arg];
      AS_UTL_loadFileList(corInputsFile, corInputs);

    } else if (AS_UTL_fileExists(argv[arg])) {
      corInputs.push_back(argv[arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (gkpName == NULL)
    err.push_back("ERROR:  no gatekeeper store (-G) supplied.\n");
  if (corName == NULL)
    err.push_back("ERROR:  no tig store (-T) supplied.\n");
  if ((corInputs.size() == 0) && (corInputsFile == NULL))
    err.push_back("ERROR:  no input tigs supplied on command line and no -L file supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G <gkpStore> -C <corStore> [input.cns]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <gkpStore>         Path to the gatekeeper store\n");
    fprintf(stderr, "  -C <corStore>         Path to the corStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L <file-of-files>    Load the tig(s) from files listed in 'file-of-files'\n");
    fprintf(stderr, "                        (WARNING: program will succeed if this file is empty)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Load the output of falconsense into the corStore and gkpStore.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  gkStore     *gkpStore = gkStore::gkStore_open(gkpName, gkStore_extend);
  gkReadData  *readData = new gkReadData;
  tgStore     *corStore = new tgStore(corName, corVers, tgStoreModify);
  tgTig       *tig      = new tgTig;

  for (uint32 ff=0; ff<corInputs.size(); ff++) {
    errno = 0;
    FILE *TI = fopen(corInputs[ff], "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", corInputs[ff], strerror(errno)), exit(1);

    fprintf(stderr, "Reading layouts from '%s'.\n", corInputs[ff]);

    while (tig->loadFromStreamOrLayout(TI) == true) {
      uint32  rID  = tig->tigID();
      gkRead *read = gkpStore->gkStore_getRead(rID);

      if (tig->consensusExists() == false)                               //  If no consensus, don't load.
        continue;

      fprintf(stderr, "Read %u raw length %u corrected length %u\n",
              rID, read->gkRead_sequenceLength(), tig->length());

      //corStore->insertTig(tig, false);                                   //  Load the data into corStore.

      //  Load the data into gkpStore.

      gkpStore->gkStore_loadReadData(tig->tigID(), readData);            //  Load old data.
      readData->gkReadData_setBasesQuals(tig->bases(), tig->quals());    //  Insert new data.
      gkpStore->gkStore_stashReadData(readData);                         //  Write.
    }

    AS_UTL_closeFile(TI, corInputs[ff]);

    fprintf(stderr, "Reading layouts from '%s' completed.\n", corInputs[ff]);
  }

  delete tig;
  delete corStore;

  delete readData;

  gkpStore->gkStore_close();

  exit(0);
}
