
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

#ifndef UTGCNS_H
#define UTGCNS_H

class cnsParameters {
public:
  void                    closeAndCleanup(void) {
    delete seqStore;   seqStore = NULL;
    delete tigStore;   tigStore = NULL;

    if (seqReads)
      for (auto it=seqReads->begin(); it != seqReads->end(); it++)
        delete it->second;

    delete seqReads;

    merylutil::closeFile(outResultsFile, outResultsName);
    merylutil::closeFile(outLayoutsFile, outLayoutsName);

    merylutil::closeFile(outSeqFileA, outSeqNameA);
    merylutil::closeFile(outSeqFileQ, outSeqNameQ);
  };

  char                   *seqName = nullptr;
  char                   *seqFile = nullptr;

  char                   *tigName = nullptr;
  uint32                  tigVers = UINT32_MAX;
  uint32                  tigPart = 0;

  uint32                  tigBgn  = 0;
  uint32                  tigEnd  = UINT32_MAX;

  char                   *outResultsName = nullptr;
  char                   *outLayoutsName = nullptr;
  char                   *outSeqNameA    = nullptr;
  char                   *outSeqNameQ    = nullptr;

  char                   *exportName     = nullptr;
  char                   *importName     = nullptr;
  bool                    dumpImport     = false;

  char                    algorithm = 'P';
  char                    aligner   = 'E';

  bool                    createPartitions = false;
  double                  partitionSize    = 1.00;   //  Size partitions to be 100% of the largest tig.
  double                  partitionScaling = 1.00;   //  Estimated tig length is 100% of actual tig length.
  double                  partitionReads   = 0.05;   //  5% of all reads can end up in a single partition.

  double                  errorRate      = 0.12;
  double                  errorRateMax   = 0.40;
  uint32                  errorRateMaxID = 0;

  uint32                  minOverlap   = 500;
  uint32                  minCoverage  = 0;

  uint32                  numFailures = 0;

  bool                    showResult = false;

  double                  maxCov = 0.0;
  uint32                  minLen = 0;
  uint32                  maxLen = UINT32_MAX;

  bool                    onlyUnassem = false;
  bool                    onlyContig  = false;

  bool                    noBubble    = false;
  bool                    noRepeat    = false;
  bool                    noSingleton = false;

  uint32                  verbosity = 0;

  sqStore                *seqStore = nullptr;
  std::map<uint32, sqRead *>  *seqReads = nullptr;
  tgStore                *tigStore = nullptr;

  FILE                   *outResultsFile = nullptr;
  FILE                   *outLayoutsFile = nullptr;
  FILE                   *outSeqFileA    = nullptr;
  FILE                   *outSeqFileQ    = nullptr;
};

