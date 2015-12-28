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

void SimpleAligner::align(dagcon::Alignment &aln) {
  NDalignment::NDalignResult ndaln;
  bool aligned = NDalignment::align(aln.qstr.c_str(), aln.qstr.size(), aln.tstr.c_str(), aln.tstr.size(), 100, true, ndaln);

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
