
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

#ifndef UNITIGPARTITION_H
#define UNITIGPARTITION_H

#include "runtime.H"
#include "system.H"
#include "strings.H"

#include "sqStore.H"
#include "tgStore.H"

#include <vector>

struct tigInfo {
  uint32   tigID       = 0;
  uint64   tigLength   = 0;
  uint64   tigChildren = 0;

  uint64   consensusArea   = 0;
  uint64   consensusMemory = 0;

  uint32   partition = uint32max;
};


class tigPartitioning {
public:
  void     loadTigInfo(tgStore             *tigStore);
  void     loadTigInfo(std::vector<tgTig *> &tigList);

  void     greedilyPartition(double   partitionSizeScale,
                             double   tigLengthScale,
                             double   maxReadsPer,
                             bool     verbose);

  void     outputPartitions(sqStore *seqStore, tgStore *tigStore, char const *storeName);
  void     reportPartitioning(FILE *partFile);


public:
  //  Inputs.
  uint32                     _nTigs    = 0;
  uint32                     _nReads   = 0;

  //  Data.
  std::vector<tgTig *> const _tigList;
  std::vector<tigInfo>       _tigInfo;

  uint32                     _nPartitions = 0;
};

#endif  //  UNITIG_PARTITION_H
