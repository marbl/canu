
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

#include "system.H"
#include "strings.H"

#include "kmers.H"

using namespace merylutil::kmers::v1;

#include "sqStore.H"
#include "tgStore.H"

#include "unitigConsensus.H"
#include "unitigPartition.H"

#include <set>
#include <map>

//
//  Essentially global structure containing command line parameters and a
//  few utility functions that clutter up processing algorithms.
//

class cnsParameters {
public:
  void      loadPartitionedReads(void);
  void      loadProcessList(void);
  tgTig    *copyTig(uint32 ti);
  bool      skipTig(tgTig *tig);
  void      unloadReads(void);
  void      closeAndCleanup(void);

public:
  char         *seqName = nullptr;
  char         *seqFile = nullptr;

  char         *tigName = nullptr;
  uint32        tigVers = UINT32_MAX;
  uint32        tigPart = 0;

  typedef std::set<uint32> intSet;

  intSet        processList;

  uint32        tigBgn  = 1;
  uint32        tigCur  = 0;
  uint32        tigEnd  = UINT32_MAX;

  char         *outResultsName = nullptr;
  char         *outLayoutsName = nullptr;
  char         *outSeqNameA    = nullptr;
  char         *outSeqNameQ    = nullptr;
  char         *outBAMName     = nullptr;

  char         *exportName     = nullptr;

  //ar         *markersName    = nullptr;

  char         *importName     = nullptr;
  readBuffer   *importFile     = nullptr;

  bool          dumpImport     = false;
  FILE         *dumpedLayouts  = nullptr;
  FILE         *dumpedReads    = nullptr;

  char          algorithm = 'P';
  char          aligner   = 'E';

  bool          createPartitions = false;
  double        partitionSize    = 1.00;   //  Size partitions to be 100% of the largest tig.
  double        partitionScaling = 1.00;   //  Estimated tig length is 100% of actual tig length.
  double        partitionReads   = 0.05;   //  5% of all reads can end up in a single partition.

  double        errorRate      = 0.12;
  double        errorRateMax   = 0.40;
  uint32        errorRateMaxID = 0;

  uint32        minOverlap   = 500;
  uint32        minCoverage  = 0;
  uint32        numIterations = 1;

  uint32        numFailures = 0;

  bool          showResult = false;

  double        maxCov = 0.0;
  uint32        minLen = 0;
  uint32        maxLen = UINT32_MAX;

  bool          onlyUnassem = false;
  bool          onlyContig  = false;

  bool          noBubble    = false;
  bool          noRepeat    = false;
  bool          noSingleton = false;

  uint32        verbosity = 0;

  typedef std::map<uint32, sqRead *> readMap;

  sqStore      *seqStore = nullptr;
  readMap       seqReads;
  tgStore      *tigStore = nullptr;

  //uint32            markersPeak = 0;
  //merylExactLookup *markers     = nullptr;

  FILE         *outResultsFile = nullptr;
  FILE         *outLayoutsFile = nullptr;
  FILE         *outSeqFileA    = nullptr;
  FILE         *outSeqFileQ    = nullptr;
  sam_hdr_t    *outBAMhp       = nullptr;
  samFile      *outBAMfp       = nullptr;
};

void   createPartitions(cnsParameters  &params);
void   exportTigs(cnsParameters  &params);

void   processTigs(cnsParameters  &params);

#endif  //  UTGCNS_H
