
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
processImportedTigs(cnsParameters  &params) {
  readBuffer *importFile      = new readBuffer(params.importName);
  FILE       *importedLayouts = nullptr;
  FILE       *importedReads   = nullptr;

  if (params.dumpImport) {
    importedLayouts = merylutil::openOutputFile(params.importName, '.', "layout");
    importedReads   = merylutil::openOutputFile(params.importName, '.', "fasta");
  }

  tgTig                       *tig = new tgTig();
  std::map<uint32, sqRead *>   reads;

  sqRead_defaultVersion = sqRead_raw | sqRead_normal;

  while (tig->importData(importFile, reads, importedLayouts, importedReads) == true) {
    if ((params.tigBgn <= tig->tigID()) &&
        (tig->tigID()  <= params.tigEnd)) {
      fprintf(stdout, "%7u %9u %7u%s", tig->tigID(), tig->length(), tig->numberOfChildren(),
              (tig->numberOfChildren() == 1) ? "\n" : "");

      //  Stash excess coverage.  Singletons report no logging.

      tgTigStashed S;

      tig->stashContains(params.maxCov, S);

      if (S.nBack > 0)
        fprintf(stdout, "  %8u %7.2fx %8u %7.2fx  %8u %7.2fx\n",
                S.nCont, (double)S.bCont / tig->length(),
                S.nStsh, (double)S.bStsh / tig->length(),
                S.nBack, (double)S.bBack / tig->length());

      //  Compute!

      tig->_utgcns_verboseLevel = params.verbosity;

      unitigConsensus  *utgcns  = new unitigConsensus(params.seqStore, params.errorRate, params.errorRateMax, params.errorRateMaxID, params.minOverlap, params.minCoverage);
      bool              success = utgcns->generate(tig, params.algorithm, params.aligner, &reads);

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

      delete utgcns;
      utgcns = nullptr;
    }
    else {
      fprintf(stdout, "%7u %9u %7u  skipped.\n", tig->tigID(), tig->length(), tig->numberOfChildren());
    }

    //  Tidy up for the next tig.

    for (auto it=reads.begin(); it != reads.end(); ++it) {
      delete it->second;
    }
    reads.clear();

    delete tig;
    tig = new tgTig();    //  Next loop needs an existing empty layout.
  }

  delete tig;

  merylutil::closeFile(importedReads);
  merylutil::closeFile(importedLayouts);

  delete importFile;
}
