
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
 *    Sergey Koren beginning on 2015-DEC-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <vector>
#include <stdint.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "Alignment.H"
#include "SimpleAligner.H"
#include "assert.h"

SimpleAligner::SimpleAligner() {
}

void SimpleAligner::align(dagcon::Alignment &aln, double errorRate) {
  NDalignment::NDalignResult ndaln;
  bool aligned = NDalignment::align(aln.qstr.c_str(), aln.qstr.size(), aln.tstr.c_str(), aln.tstr.size(), 150, true, ndaln);

  if (((double) ndaln._dist / (double) ndaln._size) > errorRate) {
     aligned = false;
  }

  if (aligned) {
     aln.start += ndaln._tgt_bgn;
     aln.end = aln.start + ndaln._tgt_end;
     aln.start++;
     aln.qstr = std::string(ndaln._qry_aln_str);
     aln.tstr = std::string(ndaln._tgt_aln_str);
  } else {
     aln.start = aln.end = 0;
     aln.qstr = std::string();
     aln.tstr = std::string();
  }

  assert(aln.qstr.length() == aln.tstr.length());
}
