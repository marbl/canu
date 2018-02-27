
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

#include "gkStore.H"

#include "clearRangeFile.H"



int
main (int argc, char **argv) {
  char            *gkpName       = NULL;
  char            *clrName       = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      clrName = argv[++arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (gkpName == NULL)
    err.push_back("ERROR:  no gatekeeper store (-G) supplied.\n");
  if (clrName == NULL)
    fprintf(stderr, "Warning:  no clear range file (-c) supplied, using full read length.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G <gkpStore> -c <clearRangeFile>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <gkpStore>         Path to the gatekeeper store\n");
    fprintf(stderr, "  -c <clearRangeFile>   Path to the file of clear ranges\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Loads results of read trimming into gkpStore.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  gkStore        *gkpStore = gkStore::gkStore_open(gkpName, gkStore_extend);
  uint32          numReads  = gkpStore->gkStore_getNumReads();
  uint32          numLibs   = gkpStore->gkStore_getNumLibraries();

  clearRangeFile *clrRange = NULL;
  if (clrName != NULL)
     clrRange = new clearRangeFile(clrName, gkpStore);


  for (uint32 rid=1; rid<=numReads; rid++)
    if (clrRange != NULL) {
       if (clrRange->isDeleted(rid) == false) {
          gkpStore->gkStore_setClearRange(rid, clrRange->bgn(rid), clrRange->end(rid));
       }
    } else {
       gkRead* read = gkpStore->gkStore_getRead(rid);
       gkpStore->gkStore_setClearRange(rid, 0, read->gkRead_sequenceLength());
    }

  delete clrRange;

  gkpStore->gkStore_close();

  exit(0);
}
