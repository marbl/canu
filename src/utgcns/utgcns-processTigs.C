
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

#include "utgcns.H"



void
processTigs(cnsParameters  &params) {
  uint32   nTigs       = 0;
  uint32   nSingletons = 0;
  uint32   numFailures = 0;


void
printHeader(cnsParameters  &params) {
  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Computing consensus for b=" F_U32 " to e=" F_U32 " with errorRate %0.4f (max %0.4f) and minimum overlap " F_U32 "\n",
          params.tigBgn, params.tigEnd, params.errorRate, params.errorRateMax, params.minOverlap);
  fprintf(stderr, "--\n");
  fprintf(stdout, "                           ----------CONTAINED READS----------  -DOVETAIL  READS-\n");
  fprintf(stdout, "  tigID    length   reads      used coverage  ignored coverage      used coverage\n");
  fprintf(stdout, "------- --------- -------  -------- -------- -------- --------  -------- --------\n");
}


  //  Load the partition file, if it exists.

  std::set<uint32>   processList = loadProcessList(params.tigName, params.tigPart);

  //  Load the partitioned reads, if they exist.

  params.seqReads = loadPartitionedReads(params.seqFile);

  //  Loop over all tigs, loading each one and processing if requested.

  for (uint32 ti=params.tigBgn; ti<=params.tigEnd; ti++) {

    if ((processList.size() > 0) &&       //  Ignore tigs not in our partition.
        (processList.count(ti) == 0))     //  (if a partition exists)
      continue;

    tgTig *tig = params.tigStore->loadTig(ti);

    if ((tig == NULL) ||                  //  Ignore non-existent and
        (tig->numberOfChildren() == 0))   //  empty tigs.
      continue;

    //  Skip stuff we want to skip.

    if (((params.onlyUnassem == true) && (tig->_class != tgTig_unassembled)) ||
        ((params.onlyContig  == true) && (tig->_class != tgTig_contig)) ||
        ((params.noSingleton == true) && (tig->numberOfChildren() == 1)) ||
        (tig->length() < params.minLen) ||
        (tig->length() > params.maxLen))
      continue;

    //  Skip repeats and bubbles.

    if (((params.noRepeat == true) && (tig->_suggestRepeat == true)) ||
        ((params.noBubble == true) && (tig->_suggestBubble == true)))
      continue;

    //  Log that we're processing.

    if (tig->numberOfChildren() > 1) {
      fprintf(stdout, "%7u %9u %7u", tig->tigID(), tig->length(), tig->numberOfChildren());
    }

    //  Stash excess coverage.  Singletons report no logging.

    tgTigStashed S;

    tig->stashContains(params.maxCov, S);

    if (S.nBack > 0) {
      nTigs++;
      fprintf(stdout, "  %8u %7.2fx %8u %7.2fx  %8u %7.2fx\n",
              S.nCont, (double)S.bCont / tig->length(),
              S.nStsh, (double)S.bStsh / tig->length(),
              S.nBack, (double)S.bBack / tig->length());
    } else {
      nSingletons++;
    }

    //  Compute!

    tig->_utgcns_verboseLevel = params.verbosity;

    unitigConsensus  *utgcns  = new unitigConsensus(params.seqStore, params.errorRate, params.errorRateMax, params.errorRateMaxID, params.minOverlap, params.minCoverage);
    bool              success = utgcns->generate(tig, params.algorithm, params.aligner, params.seqReads);

    //  Show the result, if requested.

    if (params.showResult)
      tig->display(stdout, params.seqStore, 200, 3);

    //  Unstash.

    tig->unstashContains();

    //  Save the result.

    if (params.outResultsFile)   tig->saveToStream(params.outResultsFile);
    if (params.outLayoutsFile)   tig->dumpLayout(params.outLayoutsFile);
    if (params.outSeqFileA)      tig->dumpFASTA(params.outSeqFileA);
    if (params.outSeqFileQ)      tig->dumpFASTQ(params.outSeqFileQ);

    //  Count failure.

    if (success == false) {
      fprintf(stderr, "unitigConsensus()-- tig %d failed.\n", tig->tigID());
      numFailures++;
    }

    //  Tidy up for the next tig.

    delete utgcns;        //  No real reason to keep this until here.

    params.tigStore->unloadTig(tig->tigID(), true);  //  Tell the store we're done with it
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "Processed %u tig%s and %u singleton%s.\n",
          nTigs, (nTigs == 1)             ? "" : "s",
          nSingletons, (nSingletons == 1) ? "" : "s");
  fprintf(stdout, "\n");

  if (numFailures) {
    fprintf(stderr, "WARNING:  %u tig%s failed.\n", numFailures, (numFailures == 1) ? "" : "s");
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");
  } else {
    fprintf(stderr, "Consensus finished successfully.\n");
  }
}
