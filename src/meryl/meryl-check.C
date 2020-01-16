
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kmers.H"
#include "sequence.H"
#include "bits.H"




int
main(int argc, char **argv) {
  char   *inputSeqName = NULL;
  char   *inputDBname  = NULL;
  uint64  minV         = 0;
  uint64  maxV         = UINT64_MAX;
  uint32  threads      = 1;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-sequence") == 0) {   //  INPUT READS and RANGE TO PROCESS
      inputSeqName = argv[++arg];

    } else if (strcmp(argv[arg], "-mers") == 0) {
      inputDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      minV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      maxV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      threads = strtouint32(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputSeqName == NULL)
    err.push_back("No input sequences (-sequence) supplied.\n");
  if (inputDBname == NULL)
    err.push_back("No query meryl database (-mers) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }



  map<kmer,uint32>   check;

  //  Open a database, load the kmers and values into 'check'.

  fprintf(stderr, "Open meryl database '%s'.\n", inputDBname);
  kmerCountFileReader   *merylDB    = new kmerCountFileReader(inputDBname);

  fprintf(stderr, "Convert to lookup table.\n");
  //kmerCountExactLookup  *kmerLookup = new kmerCountExactLookup(merylDB, minV, maxV);

  fprintf(stderr, "Create mapping to value.\n");
  uint64                 nKmers     = 0;

  while (merylDB->nextMer() == true) {
    kmer    kmer  = merylDB->theFMer();
    uint32  value = merylDB->theValue();

    check[kmer] = value;

    nKmers++;

    if ((nKmers % 100000) == 0) {
      fprintf(stderr, "Loaded %li kmers.\n", nKmers);
    }
  }

  delete merylDB;
  //delete kmerLookup;

  fprintf(stderr,"Loaded %lu kmers into check map of size %lu\n", nKmers, check.size());

  //

  fprintf(stderr, "Stream kmers from '%s'.\n", inputSeqName);

  dnaSeqFile  *seqFile    = new dnaSeqFile(inputSeqName);

  {
  uint32   nameMax = 0;
  char    *name    = NULL;
  uint64   seqLen  = 0;
  uint64   seqMax  = 0;
  char    *seq     = NULL;
  uint8   *qlt     = NULL;

  char     fString[64];
  char     rString[64];

  while (seqFile->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen)) {
    kmerIterator  kiter(seq, seqLen);

    while (kiter.nextMer()) {
      kmer     fMer  = kiter.fmer();
      kmer     rMer  = kiter.rmer();
      uint64   value = 0;

      if (fMer < rMer)
        value = check[fMer]--;
      else
        value = check[rMer]--;

      if (value == 0)
        fprintf(stdout, "%s\t%s\t%s ZERO\n",
                name,
                kiter.fmer().toString(fString),
                kiter.rmer().toString(rString));

    }
  }

  delete [] name;
  delete [] seq;
  delete [] qlt;
  }

  delete seqFile;

  //  Check that all values are zero.

  for (map<kmer,uint32>::iterator it=check.begin(); it != check.end(); it++) {
    kmer    k = it->first;
    uint32  v = it->second;

    if (v != 0) {
      char   kmerString[64];

      fprintf(stderr, "%s\t%u\n", k.toString(kmerString), v);
    }
  }

  exit(0);
}
