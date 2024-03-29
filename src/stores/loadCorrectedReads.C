
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

#include "system.H"
#include "strings.H"

#include "sqStore.H"
#include "tgStore.H"

#include <vector>


int
main (int argc, char **argv) {
  char const                *seqName        = NULL;
  char const                *corName        = NULL;
  int32                      corVers        = 1;

  stringList                 corInputs;
  char const                *corInputsFile  = NULL;

  bool                       updateCorStore = false;

  argc = AS_configure(argc, argv, 1);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corName = argv[++arg];
      //corVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      corInputsFile = argv[++arg];
      corInputs.load(corInputsFile);

    } else if (strcmp(argv[arg], "-u") == 0) {
      updateCorStore = true;

    } else if (fileExists(argv[arg])) {
      corInputs.add(argv[arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (seqName == NULL)
    err.push_back("ERROR:  no sequence store (-S) supplied.\n");
  if (corName == NULL)
    err.push_back("ERROR:  no tig store (-T) supplied.\n");
  if ((corInputs.size() == 0) && (corInputsFile == NULL))
    err.push_back("ERROR:  no input tigs supplied on command line and no -L file supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S <seqStore> -C <corStore> [input.cns]\n", argv[0]);
    fprintf(stderr, "  Load the output of falconsense into the corStore and seqStore.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S <seqStore>         Path to a sequence store\n");
    fprintf(stderr, "  -C <corStore>         Path to a correction store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L <file-of-files>    Load the tig(s) from files listed in 'file-of-files'\n");
    fprintf(stderr, "                        (WARNING: program will succeed if this file is empty)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -u                    Also load the populated tig layout into version 2 of the corStore.\n");
    fprintf(stderr, "                          Doesn't work; falconsense isn't updating the layout with\n");
    fprintf(stderr, "                          the position of the evidence read in the corrected read.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sqStore          *seqStore = new sqStore(seqName, sqStore_extend);
  sqRead           *read     = new sqRead;
  sqReadDataWriter *rdw      = new sqReadDataWriter(NULL);
  tgStore          *corStore = new tgStore(corName, corVers, tgStoreModify);
  tgTig            *tig      = new tgTig;

  uint64            nSkip    = 0;
  uint64            nSkipTot = 0;

  uint64            nLoad    = 0;
  uint64            nLoadTot = 0;

  fprintf(stdout, "     read       raw corrected\n");
  fprintf(stdout, "       id    length    length\n");
  fprintf(stdout, "--------- --------- ---------\n");

  fprintf(stderr, "\n");
  fprintf(stderr, "   loaded   skipped                          input file\n");
  fprintf(stderr, "--------- --------- -----------------------------------\n");

  for (uint32 ff=0; ff<corInputs.size(); ff++) {
    nSkip = 0;
    nLoad = 0;

    FILE *TI = merylutil::openInputFile(corInputs[ff]);

    while (tig->loadFromStreamOrLayout(TI) == true) {
      uint32  rID  = tig->tigID();

      if (tig->consensusExists() == false) {
        nSkip++;
        continue;
      }

      nLoad++;

      //  Load the data into corStore.

      if (updateCorStore == true)
        corStore->insertTig(tig, false);

      //  Load the data into seqStore.

      seqStore->sqStore_getRead(rID, read);                      //  Load old data for the read.

      rdw->sqReadDataWriter_importData(read);                    //  Import it into the writer.
      rdw->sqReadDataWriter_setCorrectedBases(tig->bases(),
                                              tig->length());    //  Add the corrected read.

      seqStore->sqStore_addRead(rdw);                            //  Write combined data.

      //  Log it.

      fprintf(stdout, "%9u %9u %9u %9u\n",
              rID,
              seqStore->sqStore_getReadLength(rID, sqRead_raw),
              seqStore->sqStore_getReadLength(rID, sqRead_corrected),
              tig->length());

      if (seqStore->sqStore_getReadLength(rID, sqRead_corrected) != tig->length())
        fprintf(stderr, "Read length %u differs from tig length %u\n",
                seqStore->sqStore_getReadLength(rID, sqRead_corrected),
                tig->length());
      assert(seqStore->sqStore_getReadLength(rID, sqRead_corrected) == tig->length());
    }

    merylutil::closeFile(TI, corInputs[ff]);

    fprintf(stderr, "%9" F_U64P " %9" F_U64P " %35s\n", nLoad, nSkip, corInputs[ff]);

    nSkipTot += nSkip;
    nLoadTot += nLoad;
  }

  delete tig;
  delete corStore;

  delete read;

  delete seqStore;

  fprintf(stderr, "--------- --------- -----------------------------------\n");
  fprintf(stderr, "%9" F_U64P " %9" F_U64P " %35" F_U64P "\n", nLoadTot, nSkipTot, corInputs.size());
  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}
