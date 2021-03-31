
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

#ifndef TGTIGSIZEANALYSIS
#define TGTIGSIZEANALYSIS

#include "tgTig.H"

#include <vector>

class tgTigSizeAnalysis {
public:
  tgTigSizeAnalysis(uint64 genomeSize);
  ~tgTigSizeAnalysis();

  void         evaluateTig(tgTig *tig);
  void         finalize(void);

  void         printSummary(FILE *out, char const *description, std::vector<uint32> &data);
  void         printSummary(FILE *out);

private:
  uint64               genomeSize;

  std::vector<uint32>  lenSuggestRepeat;
  std::vector<uint32>  lenSuggestCircular;

  std::vector<uint32>  lenUnassembled;
  std::vector<uint32>  lenBubble;
  std::vector<uint32>  lenContig;
};


#endif  //  TGTIGSIZEANALYSIS
