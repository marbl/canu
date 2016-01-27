
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
 *  This file is derived from:
 *
 *    src/AS_CNS/MultiAlignSizeAnalysis.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-MAR-26 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-22
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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

  void         evaluateTig(tgTig *tig, bool useGapped=true);
  void         finalize(void);

  void         printSummary(FILE *out, char *description, vector<uint32> &data);
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
