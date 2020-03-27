
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

using namespace std;


class tgTigSizeAnalysis {
public:
  tgTigSizeAnalysis(uint64 genomeSize);
  ~tgTigSizeAnalysis();

  void         evaluateTig(tgTig *tig);
  void         finalize(void);

  void         printSummary(FILE *out, char const *description, vector<uint32> &data);
  void         printSummary(FILE *out);

private:
  uint64          genomeSize;

  vector<uint32>  lenSuggestRepeat;
  vector<uint32>  lenSuggestCircular;

  vector<uint32>  lenUnassembled;
  vector<uint32>  lenBubble;
  vector<uint32>  lenContig;
};


#endif  //  TGTIGSIZEANALYSIS
