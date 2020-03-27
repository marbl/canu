
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

#include "clearRangeFile.H"

#include <vector>
using namespace std;


int
main (int argc, char **argv) {
  char const            *seqName       = NULL;

  vector<char const *>   clrName;
  vector<uint32>         bgnID;
  vector<uint32>         endID;

  char const            *outName       = NULL;

  bool                   verbose       = false;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      clrName.push_back(argv[++arg]);
      bgnID.push_back(strtouint32(argv[++arg]));
      endID.push_back(strtouint32(argv[++arg]));

    } else if (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-v") == 0) {
      verbose = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqName == NULL)
    err.push_back("ERROR:  no sequence store (-S) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S <seqStore> -c <bgnID> <endID> <clearRangeFile> -o <clearRangeFile>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S <seqStore>                    Path to the sequence store\n");
    fprintf(stderr, "  -c <clearRangeFile> <bgn> <end>  Path to the file of clear ranges,\n");
    fprintf(stderr, "                                   along with the (inclusive) range of\n");
    fprintf(stderr, "                                   read IDs that have clear ranges set\n");
    fprintf(stderr, "  -o <clearRangeFile>              Path to output clear ranges.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -v                    Report clear range changes to stderr\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Merges multiple clear range files into one.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sqStore        *seqStore = new sqStore(seqName, sqStore_extend);
  uint32          numReads  = seqStore->sqStore_lastReadID();
  uint32          numLibs   = seqStore->sqStore_lastLibraryID();

  clearRangeFile *outRange = new clearRangeFile(outName, seqStore);

  for (uint32 ii=0; ii<clrName.size(); ii++) {
    clearRangeFile *clrRange = new clearRangeFile(clrName[ii], seqStore);

    for (uint32 rid=bgnID[ii]; rid<=endID[ii]; rid++) {
      if (verbose == true)
        fprintf(stderr, "%u\t%7u-%-7u\t%7u-%-7u\n",
                rid,
                seqStore->sqStore_getClearBgn(rid), seqStore->sqStore_getClearEnd(rid),
                clrRange->bgn(rid), clrRange->end(rid));

      outRange->setbgn(rid) = clrRange->bgn(rid);
      outRange->setend(rid) = clrRange->end(rid);
    }

    delete clrRange;
  }

  delete outRange;

  delete seqStore;

  exit(0);
}
