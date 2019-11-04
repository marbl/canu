
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
 *    Sergey Koren beginning on 2017-OCT-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "sqStore.H"

#include "clearRangeFile.H"



int
main (int argc, char **argv) {
  char            *seqName = NULL;
  char            *clrName = NULL;

  sqStore_mode     sqMode  = sqStore_extend;

  bool             verbose = false;
  bool             modify  = true;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      clrName = argv[++arg];

    } else if (strcmp(argv[arg], "-v") == 0) {
      verbose = true;

    } else if (strcmp(argv[arg], "-n") == 0) {
      modify = false;
      sqMode = sqStore_readOnly;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqName == NULL)
    err.push_back("ERROR:  no sequence store (-S) supplied.\n");
  if (clrName == NULL)
    fprintf(stderr, "Warning:  no clear range file (-c) supplied, using full read length.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S <seqStore> -c <clearRangeFile>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S <seqStore>         Path to the sequence store\n");
    fprintf(stderr, "  -c <clearRangeFile>   Path to the file of clear ranges\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -v                    Report clear range changes to stderr\n");
    fprintf(stderr, "  -n                    Don't apply changes\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Loads results of read trimming into seqStore.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sqStore        *seqStore = sqStore::sqStore_open(seqName, sqMode);
  uint32          numReads  = seqStore->sqStore_getNumReads();
  uint32          numLibs   = seqStore->sqStore_getNumLibraries();

  //  If no clear range file, set clear range to the entire read.

  if (clrName == NULL) {
    for (uint32 rid=1; rid<=numReads; rid++) {
      sqRead* read = seqStore->sqStore_getRead(rid);

      seqStore->sqStore_setClearRange(rid, 0, read->sqRead_sequenceLength());
    }
  }

  //  Otherwise, set to whatever the clear range file says.

  else {
    clearRangeFile *clrRange = new clearRangeFile(clrName, seqStore);

    for (uint32 rid=1; rid<=numReads; rid++) {
      sqRead* read = seqStore->sqStore_getRead(rid);

      if (verbose == true)
        fprintf(stderr, "%u\t%7u-%-7u\t%7u-%-7u\n",
                rid,
                read->sqRead_clearBgn(), read->sqRead_clearEnd(),
                clrRange->bgn(rid), clrRange->end(rid));

      if (clrRange->isDeleted(rid) == true)
        continue;

      if (modify == true)
        seqStore->sqStore_setClearRange(rid, clrRange->bgn(rid), clrRange->end(rid));
    }

    delete clrRange;
  }

  seqStore->sqStore_close();

  exit(0);
}
