
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




tgTig *
loadTigFromImport(cnsParameters &params) {
  tgTig *tig = nullptr;

 tryImportAgain:
  tig = new tgTig;
 
  params.unloadReads();                       //  Forget any reads from the last tig.

  if (tig->importData(params.importFile,      //  Load the next tig/reads from the package.
                      params.seqReads,        //  If no next, we're done.
                      params.dumpedLayouts,
                      params.dumpedReads) == false) {
    delete tig;
    return nullptr;
  }

  if (params.skipTig(tig)) {                  //  If params say to skip the tig,
    delete tig;                               //  forget the tig and
    goto tryImportAgain;                      //  load another one.
  }

  return tig;
}



tgTig *
loadTigFromStore(cnsParameters &params) {
  tgTig *tig = nullptr;

 tryLoadAgain:
  if (params.tigCur < params.tigBgn)    //  Advance to the first or next tig.  On the first
    params.tigCur = params.tigBgn;      //  call, tigCur is 0, and we'll either set it to
  else                                  //  tigBgn (above) or to 1 (below).  On later calls,
    params.tigCur++;                    //  we'll always increment tigCur.
  
  if (params.tigCur > params.tigEnd)   //  If current is past the end, we're done.
    return nullptr;

  tig = params.copyTig(params.tigCur);  //  Otherwise, load A COPY of the tig from the store.

  if (params.skipTig(tig)) {            //  If params say to skip the tig,
    delete tig;                         //  forget the tig (we own the copy) and
    goto tryLoadAgain;                  //  load another one.
  }

  return tig;
}



tgTig *
loadNextTig(cnsParameters &params) {
  tgTig  *tig = nullptr;

  if (params.importFile)
    tig = loadTigFromImport(params);
  else
    tig = loadTigFromStore(params);

  return tig;
}



void
processTigs(cnsParameters  &params) {
  uint32      nTigs           = 0;
  uint32      nSingletons     = 0;
  uint32      numFailures     = 0;

  //  Print a lovely header for the progress report.

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Computing consensus for b=" F_U32 " to e=" F_U32 " with errorRate %0.4f (max %0.4f) and minimum overlap " F_U32 "\n",
          params.tigBgn, params.tigEnd, params.errorRate, params.errorRateMax, params.minOverlap);
  fprintf(stderr, "--\n");
  fprintf(stdout, "                           ----------CONTAINED READS----------  -DOVETAIL  READS-\n");
  fprintf(stdout, "  tigID    length   reads      used coverage  ignored coverage      used coverage\n");
  fprintf(stdout, "------- --------- -------  -------- -------- -------- --------  -------- --------\n");

  //  Load the partitioned reads or open the package.

  if (params.importName) {
    params.importFile    = new readBuffer(params.importName);
    params.dumpedLayouts = merylutil::openOutputFile(params.importName, '.', "layout", params.dumpImport);
    params.dumpedReads   = merylutil::openOutputFile(params.importName, '.', "fasta",  params.dumpImport);
  }
  else {
    params.loadPartitionedReads();
  }

  //  Loop over all tigs, loading each one and processing if requested.

  for (tgTig *tig=loadNextTig(params); tig != nullptr; tig=loadNextTig(params)) {

    //  Log that we're processing (but ignore singletons) and filter
    //  contained and ignore reads.

    nTigs       += (tig->numberOfChildren() > 1) ? 1 : 0;
    nSingletons += (tig->numberOfChildren() > 1) ? 0 : 1;

    tig->filterContains(params.maxCov, false);

    if (tig->numberOfChildren() > 1)
      fprintf(stdout, "  %8lu %7.2fx %8lu %7.2fx  %8lu %7.2fx\n",  //  The start of this line
              tig->nStashCont(), tig->cStashCont(),                //  is printed by
              tig->nStashStsh(), tig->cStashStsh(),                //  cnsParameters::skipTig().
              tig->nStashBack(), tig->cStashBack());

    //  Compute!
    unitigConsensus  *utgcns  = nullptr;
    bool              success = false;

    utgcns  = new unitigConsensus(params.seqStore,
                                  params.errorRate, params.errorRateMax, params.errorRateMaxID,
                                  params.minOverlap,
                                  params.minCoverage);

    success = utgcns->generate(tig,
                               params.algorithm,
                               params.aligner,
                               params.seqReads);

    delete utgcns;

    //  Show the result, if requested.

    if (params.showResult)
      tig->display(stdout, params.seqStore, 200, 3);

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

    delete tig;
  }

  fprintf(stdout, "------- --------- -------  -------- -------- -------- --------  -------- --------\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "-- Processed %u tig%s and %u singleton%s.\n",
          nTigs, (nTigs == 1)             ? "" : "s",
          nSingletons, (nSingletons == 1) ? "" : "s");
  fprintf(stdout, "-- \n");

  if (numFailures) {
    fprintf(stderr, "-- WARNING:  %u tig%s failed.\n", numFailures, (numFailures == 1) ? "" : "s");
    fprintf(stderr, "-- \n");
    fprintf(stderr, "-- Consensus did NOT finish successfully.\n");
  } else {
    fprintf(stderr, "-- Consensus finished successfully.\n");
  }
}
