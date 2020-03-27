
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

#include "sweatShop.H"
#include "system.H"
#include "sequence.H"

#include <pthread.h>

#include "sqStore.H"
#include "sqCache.H"
#include "ovStore.H"

#include "alignStats.H"
#include "overlapAlign-globalData.H"
#include "overlapAlign-threadData.H"
#include "overlapAlign-computation.H"

#include "clearRangeFile.H"



void *
overlapReader(void *G) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = NULL;

  while ((g->curID <= g->endID) &&                    //  Skip any reads with no overlaps.
         (g->ovlStore->numOverlaps(g->curID) == 0))
    g->curID++;

  if (g->curID <= g->endID) {                         //  Make a new computation object,
    s = new maComputation(g->curID,                   //  and advance to the next read.
                          g->readData,
                          g->seqCache,
                          g->ovlStore,
                          g->verboseTrim,
                          g->verboseAlign);
    g->curID++;
  }

  return(s);
}



void
overlapWriter(void *G, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = (maComputation *)S;
  uint32            n = 0;

  //  Count the number of overlaps to output.

  for (uint64 oo=0; oo<s->_overlapsLen; oo++)
    if (s->_overlaps[oo].evalue() < AS_MAX_EVALUE)
      n++;

  //  If there are overlaps, output them.

  if (n > 0) {
    fprintf(stdout, "\n");
    fprintf(stdout, "%6u %6d %6d %s\n", s->_aID, 0, s->_readData[s->_aID].trimmedLength, s->_aRead);

    for (uint64 oo=0; oo<s->_overlapsLen; oo++) {
      if (s->_overlaps[oo].evalue() == AS_MAX_EVALUE)    //  Skip garbage.
        continue;

      g->outFile->writeOverlap(&s->_overlaps[oo]);

      fprintf(stdout, "%6u %6d %6d %s\n",
              s->_overlaps[oo].b_iid,
              (int32)s->_overlaps[oo].dat.ovl.ahg5,
              (int32)s->_overlaps[oo].dat.ovl.ahg3,
              s->_alignsB[oo]);
    }
  }

  //  Cleanup after the compute.

  delete s;
}



void
trimWriter(void *G, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = (maComputation *)S;

  //  Do nothing, just delete.

  delete s;
}



void
overlapRecompute(void *G, void *T, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maThreadData     *t = (maThreadData  *)T;
  maComputation    *s = (maComputation *)S;

  //fprintf(stderr, "Processing read %u with %u overlaps.\n", s->_aID, s->_overlapsLen);

  s->computeAlignments(g->minOverlapLength, g->maxErate);
};



void
overlapTrim(void *G, void *T, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maThreadData     *t = (maThreadData  *)T;
  maComputation    *s = (maComputation *)S;

  //fprintf(stderr, "Processing read %u with %u overlaps.\n", s->_aID, s->_overlapsLen);

  //  Trim the read.
  //
  //  Note that output is set directly in the trReadData array in trGlobalData.
  //  See the _readData member in maComputation, and overlapReader() above.
  //
  s->trimRead(g->minOverlapLength, g->maxErate);
};



void
saveClearRanges(trGlobalData *g) {

  if (g->clearRangesFileName == NULL)
    return;

  clearRangeFile *cr = new clearRangeFile(g->clearRangesFileName);

  for (uint32 ii=g->bgnID; ii<=g->endID; ii++) {
    int32  bgn = g->readData[ii].clrBgn;
    int32  end = g->readData[ii].clrEnd;

    cr->setClearRange(ii, bgn, end, (end - bgn < g->minReadLength));
  }

  cr->saveData();

  delete cr;
}



void
loadClearRanges(trGlobalData *g, vector<char const *> &trimFiles) {

  //  For each clear range file, load the data, then copy to
  //  our storage.

  for (uint32 ff=0; ff<trimFiles.size(); ff++) {
    fprintf(stderr, "Loading clear ranges from '%s'\n", trimFiles[ff]);

    clearRangeFile *cr = new clearRangeFile(trimFiles[ff]);

    cr->loadData();

    for (uint32 ii=g->bgnID; ii<=g->endID; ii++) {
      int32  bgn = cr->bgn(ii);
      int32  end = cr->end(ii);

      if (cr->valid(ii))
        g->readData[ii].setClear(bgn, end);
    }

    delete cr;
  }
}



void
setClearRanges(trGlobalData *g) {

  for (uint32 ii=1; ii<g->seqStore->sqStore_lastReadID()+1; ii++) {
    int32       bgn = g->readData[ii].clrBgn;
    int32       end = g->readData[ii].clrEnd;

    if (end - bgn < 1000) {
      fprintf(stderr, "read %6u length %6u clear %6u %6u -- DELETE\n",
              ii,
              g->seqStore->sqStore_getReadLength(ii),
              g->readData[ii].clrBgn, g->readData[ii].clrEnd);
      g->seqStore->sqStore_setIgnored(ii, false, true);
      continue;
    }

    fprintf(stderr, "read %6u length %6u clear %6u %6u\n",
            ii,
            g->seqStore->sqStore_getReadLength(ii),
            g->readData[ii].clrBgn, g->readData[ii].clrEnd);

    assert(g->readData[ii].clrBgn <  g->readData[ii].clrEnd);
    assert(g->readData[ii].clrEnd <= g->seqStore->sqStore_getReadLength(ii));

    g->seqStore->sqStore_setClearRange(ii, g->readData[ii].clrBgn, g->readData[ii].clrEnd, false);
  }

  sqRead_setDefaultVersion(sqRead_raw | sqRead_trimmed);
}



void
alignOverlaps(trGlobalData *g, bool isTrimming) {

  //  Set the range of overlaps to process.

  g->resetOverlapIteration();

  //  If only one thread, don't use sweatShop.  Easier to debug
  //  and works with valgrind.

  if (g->numThreads == 1) {
    maThreadData  *t = new maThreadData(g, 0);

    while (1) {
      maComputation *c = (maComputation *)overlapReader(g);

      if (c == NULL)
        break;

      if (isTrimming) {
        overlapTrim(g, t, c);
        trimWriter(g, c);
      }

      else {
        overlapRecompute(g, t, c);
        overlapWriter(g, c);
      }
    }

    delete t;
  }

  //  Use all the CPUs!

  else {
    maThreadData **td = new maThreadData * [g->numThreads];
    sweatShop     *ss = NULL;

    if (isTrimming) {
      ss = new sweatShop(overlapReader, overlapTrim, trimWriter);
    }

    else {
      ss = new sweatShop(overlapReader, overlapRecompute, overlapWriter);
    }

    ss->setLoaderQueueSize(512);
    ss->setWriterQueueSize(16 * 1024);    //  Otherwise skipped reads hold up the queue.

    ss->setNumberOfWorkers(g->numThreads);

    for (uint32 w=0; w<g->numThreads; w++)
      ss->setThreadData(w, td[w] = new maThreadData(g, w));

    //  Turn on progress reports if debugging output is disabled.
    ss->run(g, (isTrimming == true) ? (g->verboseTrim == 0) : (g->verboseAlign == 0));

    delete ss;

    for (uint32 w=0; w<g->numThreads; w++)
      delete td[w];

    delete [] td;
  }
}



int
main(int argc, char **argv) {
  trGlobalData   *g = new trGlobalData;

  vector<char const *>  trimFiles;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0)
      g->seqStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      g->ovlStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-r") == 0)
      decodeRange(argv[++arg], g->bgnID, g->endID);



    else if (strcmp(argv[arg], "-threads") == 0)
      g->numThreads = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-memory") == 0)
      g->memLimit = atoi(argv[++arg]);



    else if (strcmp(argv[arg], "-overlap-erate") == 0)
      g->maxErate = atof(argv[++arg]);

    else if (strcmp(argv[arg], "-overlap-length") == 0)
      g->minOverlapLength = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-read-length") == 0)
      g->minReadLength = atoi(argv[++arg]);



    else if (strcmp(argv[arg], "-trim") == 0) {
      g->clearRangesFileName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-apply") == 0) {
      while ((arg+1 < argc) &&
             (argv[arg+1][0] != '-'))
        trimFiles.push_back(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-align") == 0) {
      g->outFileName = argv[++arg];
    }



    else if (strcmp(argv[arg], "-V") == 0) {
      g->verboseTrim++;
      g->verboseAlign++;
    }

    else if (strcmp(argv[arg], "-Vt") == 0)
      g->verboseTrim++;

    else if (strcmp(argv[arg], "-Va") == 0)
      g->verboseAlign++;

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

#if 0
  if (g->seqStoreName == NULL)   err.push_back("No sequence store (-S option) supplied.\n");
  if (g->ovlStoreName == NULL)   err.push_back("No overlap store (-O option) supplied.\n");

  if ((g->outFileName         == NULL) &&
      (g->clearRangesFileName == NULL))
    err.push_back("At least one of -trim and -align must be supplied.\n");
#endif

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -S seqStore       Mandatory, path to seqStore sequences.\n");
    fprintf(stderr, "  -O ovlStore       Mandatorym path to ovlStore overlaps.\n");
    fprintf(stderr, "  -r bgnID[-endID]  Process reads bgnID to endID, inclusive.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute trimming for a subset of reads.\n");
    fprintf(stderr, "  -trim <outputName>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Apply trimming to the store.\n");
    fprintf(stderr, "  -apply <trimFile> <...>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute alignments for a subset of reads.  All reads must be trimmed prior.\n");
    fprintf(stderr, "  -trim <inputName>\n");
    fprintf(stderr, "  -align <outputName.ovlStore>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute trimming and alignments for all reads.\n");
    fprintf(stderr, "  -align <outputName.ovlStore>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Parameters:\n");
    fprintf(stderr, "  -erate e          Overlaps are computed at 'e' fraction error; must be larger than the original erate\n");
    fprintf(stderr, "  -partial          Overlaps are 'overlapInCore -S' partial overlaps\n");
    fprintf(stderr, "  -memory m         Use up to 'm' GB of memory\n");
    fprintf(stderr, "  -threads n        Use up to 'n' cores\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Advanced options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -V, -Vt, -Va      Increase debug verbosity. -Vt increases only trimming; -Va increases only alignment.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //
  //  If there is a clear range file name:
  //    and no output overlaps, compute just trimming.
  //    and    output overlaps, load the trim range, then overlap.
  //
  //  If there is not a clear range file name:
  //    and no output overlaps, error.
  //    and    output overlaps, compute trimmin then overlaps (for small genomes).
  //

  if       (trimFiles.size() > 0) {
    g->initialize(sqStore_extend);

    fprintf(stderr, "LOAD TRIMMING.\n");
    loadClearRanges(g, trimFiles);

    setClearRanges(g);
  }

  else if ((g->clearRangesFileName != NULL) &&
           (g->outFileName         == NULL)) {
    g->initialize();

    fprintf(stderr, "COMPUTE TRIMMING and STOP.\n");
    alignOverlaps(g, true);

    saveClearRanges(g);
  }

  else if ((g->clearRangesFileName == NULL) &&
           (g->outFileName         != NULL)) {
    g->initialize();

    fprintf(stderr, "ALIGNING OVERLAPS.\n");
    alignOverlaps(g, false);
  }

  else if ((g->clearRangesFileName != NULL) &&
           (g->outFileName         != NULL)) {
    g->initialize(sqStore_extend);

    fprintf(stderr, "COMPUTE TRIMMING.\n");
    alignOverlaps(g, true);

    saveClearRanges(g);
    setClearRanges(g);

    fprintf(stderr, "ALIGNING OVERLAPS.\n");
    alignOverlaps(g, false);
  }

  //  All done!

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}
