
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



void
doTest(uint32 bgn, uint32 end, char *bases, bool norm) {

  uint32   nlen = strlen(bases);

  char     *bcomp = new char   [ nlen + 1 ];
  uint32   *ntoc  = new uint32 [ nlen + 1 ];

  uint32    clen  = homopolyCompress(bases, nlen, bcomp, ntoc);

  uint32    nbgn, nend;
  uint32    cbgn, cend;

  assert(bgn < end);

  if (norm) {
    assert(end <= nlen);

    nbgn = bgn;
    nend = end;

    cbgn = ntoc[nbgn];
    cend = ntoc[nend];
  }

  else {
    assert(end <= clen);

    cbgn = bgn;
    cend = end;

    nbgn = 0;
    while (ntoc[nbgn] <= cbgn)
      nbgn++;
    nbgn--;

    nend = nlen;
    while (cend <= ntoc[nend])
      nend--;
    nend++;
  }

  //  Probably need real bounds checking.

  assert(nbgn <= nend);
  assert(nend <= nlen);

  assert(cbgn <= cend);
  assert(cend <= clen);

  fprintf(stderr, "  normal   compressed\n");
  fprintf(stderr, "--------   ----------\n");

  for (uint32 ii=0; ii<=nlen; ii++) {
    uint32 jj = ntoc[ii];

    fprintf(stderr, "%c %4u %c %s %c %4u %c\n",
            ((ii == nbgn) || (ii == nend)) ? '-' : ' ',
            ii,
            (bases[ii] != 0) ? bases[ii] : '.',
            (bases[ii] == bcomp[jj]) ? "==" : "!!",
            (bcomp[jj] != 0) ? bcomp[jj] : '.',
            jj,
            ((jj == cbgn) || (jj == cend)) ? '-' : ' ');
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "bgn %4u %4u\n", nbgn, cbgn);
  fprintf(stderr, "end %4u %4u\n", nend, cend);

  delete [] ntoc;
  delete [] bcomp;
}



int
main (int argc, char **argv) {
  char const      *seqName = NULL;
  char const      *clrName = NULL;

  sqStore_mode     sqMode  = sqStore_extend;

  bool             verbose = false;
  bool             modify  = true;

  bool             testnorm  = false;
  bool             testcomp  = false;
  char            *testbases = NULL;
  uint32           testbgn   = 0;
  uint32           testend   = 0;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
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

    } else if (strcmp(argv[arg], "-testnorm") == 0) {
      testnorm  = true;
      testbgn   = strtouint32(argv[++arg]);
      testend   = strtouint32(argv[++arg]);
      testbases =             argv[++arg];

    } else if (strcmp(argv[arg], "-testcomp") == 0) {
      testcomp  = true;
      testbgn   = strtouint32(argv[++arg]);
      testend   = strtouint32(argv[++arg]);
      testbases =             argv[++arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((testnorm) || (testcomp)) {
    if (testbgn < testend) {
      doTest(testbgn, testend, testbases, testnorm);
    } else {
      fprintf(stderr, "Invalid clear range %u-%d\n", testbgn, testend);
    }
    exit(0);
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
    fprintf(stderr, "  -testnorm b e s       Test translating trim points between\n");
    fprintf(stderr, "  -testcomp b e s       normal and compressed sequences.  's' must\n");
    fprintf(stderr, "                        be normal (uncompressed) sequence.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Loads results of read trimming into seqStore.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sqStore        *seqStore = new sqStore(seqName, sqMode);
  uint32          numReads = seqStore->sqStore_lastReadID();
  uint32          numLibs  = seqStore->sqStore_lastLibraryID();

  //  Reset the default version to the non-trimmed version of the read.

  sqRead_setDefaultVersion(sqRead_defaultVersion & ~sqRead_trimmed);

  //  If no clear range file, set clear range to the entire read.

  if (verbose == true) {
    fprintf(stderr, "   readID   readLen  clearBgn  clearEnd    status\n");
    fprintf(stderr, "--------- --------- --------- --------- ---------\n");
  }


  if (clrName == NULL) {
    for (uint32 rid=1; rid<=numReads; rid++) {
      sqReadSeq  *rseq = seqStore->sqStore_getReadSeq(rid);

      if (rseq == NULL)
        continue;

      rseq->sqReadSeq_setClearRange(0, seqStore->sqStore_getReadLength(rid));
    }
  }

  //  Otherwise, set to whatever the clear range file says.

  else {
    clearRangeFile *clrRange = new clearRangeFile(clrName, seqStore);

    for (uint32 rid=1; rid<=numReads; rid++) {
      uint32      rlen  = seqStore->sqStore_getReadLength(rid);
      uint32      nbgn  = clrRange->bgn(rid);
      uint32      nend  = clrRange->end(rid);
      bool        bogus = false;

      //  If not a valid read, do nothing.

      if (seqStore->sqStore_isValidRead(rid) == false)
        continue;

      //  If a bogus clear range, reset to 0,0 (for logging) and ensure it is
      //  flagged for deletion.  Overlap based trimming is using UINT32_MAX
      //  as a sentinel to say 'deleted', which is gross.

      if ((nbgn >  rlen) ||
          (nend >  rlen) ||
          (nbgn >= nend)) {
        assert(clrRange->isDeleted(rid) == true);

        bogus = true;
        nbgn  = 0;
        nend  = 0;
      }

      //  If already flagged as deleted, keep the clear ranges (also for
      //  logging).

      if (clrRange->isDeleted(rid) == true)
        bogus = true;

      //  Set clear ranges for both normal and compressed versions of the read.
      //  But if bogus, flag both normal and compressed versions
      //  as ignored.

      if (modify == true)
        seqStore->sqStore_setClearRange(rid, nbgn, nend, bogus);

      //  Log the new clear range.

      if (verbose == true)
        fprintf(stderr, "%9u %9u %9u %9u%s\n",
                rid,
                rlen,
                nbgn, nend,
                bogus ? "   deleted" : "");
    }

    delete clrRange;
  }

  delete seqStore;

  exit(0);
}
