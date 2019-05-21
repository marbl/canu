
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
 *    Brian P. Walenz beginning on 2019-APR-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

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



void *
overlapReader(void *G) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = NULL;

  if (g->ovlStore) {
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
  }

  //if (g->ovlFile) {
  //  s = new maComputation(g->ovlFile);
  //}

  return(s);
}



void
overlapWriter(void *G, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = (maComputation *)S;

  for (uint64 oo=0; oo<s->_overlapsLen; oo++) {
    ovOverlap  *ovl = s->_overlaps + oo;

    if (ovl->evalue() == AS_MAX_EVALUE)
      continue;

    if (g->outStore)
      g->outStore->writeOverlap(ovl);

    if (g->outFile)
      g->outFile->writeOverlap(ovl);
  }

  if (s->_alignsLen > 0) {
    fprintf(stdout, "\n");
    fprintf(stdout, "%6u %6d %6d %s\n", s->_aID, 0, s->_readData[s->_aID].trimmedLength, s->_aRead);

    for (uint32 oo=0; oo<s->_alignsLen; oo++) {
      //fprintf(stdout, "\n");
      //fprintf(stdout, "%6u %6d %6d\n", s->_alignsOvl[oo]->b_iid, (int32)s->_alignsOvl[oo]->dat.ovl.ahg5, (int32)s->_alignsOvl[oo]->dat.ovl.ahg3);
      //fprintf(stdout, "%s\n", s->_alignsA[oo]);
      //fprintf(stdout, "%s\n", s->_alignsB[oo]);
      fprintf(stdout, "%6u %6d %6d %s\n",
              s->_alignsOvl[oo]->b_iid,
              (int32)s->_alignsOvl[oo]->dat.ovl.ahg5,
              (int32)s->_alignsOvl[oo]->dat.ovl.ahg3,
              s->_alignsB[oo]);
    }
  }

  delete s;
}



void
trimWriter(void *G, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = (maComputation *)S;

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

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0)
      g->seqStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      g->ovlStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-o") == 0)
      g->outStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-r") == 0)
      decodeRange(argv[++arg], g->bgnID, g->endID);

    else if (strcmp(argv[arg], "-t") == 0)
      g->numThreads = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-erate") == 0)
      g->maxErate = atof(argv[++arg]);

    else if (strcmp(argv[arg], "-memory") == 0)
      g->memLimit = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-len") == 0)
      g->minOverlapLength = atoi(argv[++arg]);

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

  if (g->seqStoreName == NULL)   err.push_back("No sequence store (-S option) supplied.\n");
  if (g->ovlStoreName == NULL)   err.push_back("No overlap store (-O option) supplied.\n");
  if (g->outStoreName == NULL)   err.push_back("No output store (-o option) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -S seqStore       Mandatory, path to seqStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Inputs can come from either a store or a file.\n");
    fprintf(stderr, "  -O ovlStore       \n");
    fprintf(stderr, "  -O ovlFile        \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If from an ovlStore, the range of reads processed can be restricted.\n");
    fprintf(stderr, "  -r bgnID[-endID]  process reads bgnID to endID, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs will be written to a store or file, depending on the input type\n");
    fprintf(stderr, "  -o ovlStore       \n");
    fprintf(stderr, "  -o ovlFile        \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -erate e          Overlaps are computed at 'e' fraction error; must be larger than the original erate\n");
    fprintf(stderr, "  -partial          Overlaps are 'overlapInCore -S' partial overlaps\n");
    fprintf(stderr, "  -memory m         Use up to 'm' GB of memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t n              Use up to 'n' cores\n");
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

  g->initialize();



  if (fileExists("clear-range")) {
    fprintf(stderr, "LOADING CLEAR RANGES.\n");

    FILE *clr = AS_UTL_openInputFile("clear-range");

    for (uint32 ii=0; ii<g->seqStore->sqStore_getNumRawReads()+1; ii++) {
      loadFromFile(g->readData[ii].clrBgn, "", clr);
      loadFromFile(g->readData[ii].clrEnd, "", clr);

      g->readData[ii].trimmedLength = g->readData[ii].clrEnd - g->readData[ii].clrBgn;
    }

    AS_UTL_closeFile(clr);
  }

  else {
    fprintf(stderr, "TRIMMING READS.\n");
    alignOverlaps(g, true);

    FILE *clr = AS_UTL_openOutputFile("clear-range");

    for (uint32 ii=0; ii<g->seqStore->sqStore_getNumRawReads()+1; ii++) {
      writeToFile(g->readData[ii].clrBgn, "", clr);
      writeToFile(g->readData[ii].clrEnd, "", clr);
    }

    AS_UTL_closeFile(clr);

    fprintf(stderr, "DONE.\n");
    exit(0);
  }



  //  Set all the clear ranges
  for (uint32 ii=0; ii<g->seqStore->sqStore_getNumRawReads()+1; ii++)
    g->seqStore->sqStore_setClearRange(ii, g->readData[ii].clrBgn, g->readData[ii].clrEnd);

  sqRead_setDefaultVersion(sqRead_trimmed);



  fprintf(stderr, "ALIGNING OVERLAPS.\n");
  alignOverlaps(g, false);

  //  All done!

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}
