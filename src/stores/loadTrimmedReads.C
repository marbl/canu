
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
      sqReadSeq  *rseq = seqStore->sqStore_getReadSeq(rid);
      uint32      rlen = seqStore->sqStore_getReadLength(rid);

      uint32      nbgn  = clrRange->bgn(rid);
      uint32      nend  = clrRange->end(rid);
      uint32      nlen  = nend - nbgn;

      if ((rseq == NULL) ||
          (seqStore->sqStore_isValidRead(rid) == false))
        continue;

      //  If a bogus clear range, reset to 0,0 and ensure it is flagged for
      //  deletion.  Overlap based trimming is using UINT32_MAX as a sentinel
      //  to say 'deleted', which is gross.

      if (modify == true) {
        if ((nbgn <= rlen) &&
            (nend <= rlen) &&
            (nbgn <= nend))
          rseq->sqReadSeq_setClearRange(nbgn, nend);
        else {
          nbgn = 0;  //  For display.
          nend = 0;
          assert(clrRange->isDeleted(rid) == true);
        }

        if (clrRange->isDeleted(rid) == true)
          rseq->sqReadSeq_setIgnoreT();
      }

      if (verbose == true) {
        fprintf(stderr, "%9u %9u %9u %9u%s\n",
                rid,
                rlen,
                nbgn, nend,
                clrRange->isDeleted(rid) ? "   deleted" : "");
      }
    }

    delete clrRange;
  }

  delete seqStore;

  exit(0);
}
