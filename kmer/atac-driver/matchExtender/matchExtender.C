// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <vector>
#include <algorithm>
using namespace std;

#include "bio++.H"
#include "atac.H"
#include "match.H"


bool
trim_to_pct(vector<match_s *>& matches, u32bit midx, double pct);


void
extend_match_backward(vector<match_s *>& matches,
                      u32bit midx,
		      u32bit min_start_pos);

bool
can_reach_nearby_match(match_s *src, match_s *dest);

bool
extend_match_forward(vector<match_s *>& matches, u32bit midx, match_s *target);

u32bit
extend_matches_on_diagonal(vector<match_s *>& matches, u32bit diag_start);




class MatchCompare {
public:
  int operator()(const match_s *m1, const match_s *m2) {
    return(*m1 < *m2);
  }
};





//  Read matches until the iid differs.  Leave the next match in inLine.
//
void
readMatches(char *inLine,
            FastACache  *C1,
            FastACache  *C2,
            vector<match_s *> &fwdMatches,
            vector<match_s *> &revMatches) {

  FastASequenceInCore *S1 = 0L;
  FastASequenceInCore *S2 = 0L;
  u32bit thisiid1 = ~u32bitZERO;
  u32bit thisiid2 = ~u32bitZERO;
  bool   cont;

  do {
    //chomp(inLine); fprintf(stderr, "inline=%s\n", inLine);

    u32bit  iid1=0, pos1=0, len1=0, ori1=0;
    u32bit  iid2=0, pos2=0, len2=0, ori2=0;

    cont = false;

    splitToWords *C = new splitToWords(inLine);

    if (((*C)[0][0] == 'M') && ((*C)[1][0] == 'x')) {
      decodeMatch(*C,
                  iid1, pos1, len1, ori1,
                  iid2, pos2, len2, ori2);

      //  If this is the first time through this loop (this call
      //  anyway) remember the first iid's we see.
      //
      if ((thisiid1 == ~u32bitZERO) && (thisiid2 == ~u32bitZERO)) {
        thisiid1 = iid1;
        thisiid2 = iid2;
        S1       = C1->getSequence(iid1);
        S2       = C2->getSequence(iid2);
      }

      //  If the iid's are the same as the first time, keep going, and
      //  add the match to the list.
      //
      if ((thisiid1 == iid1) && (thisiid2 == iid2)) {
        cont = true;

        if (ori1 == ori2)
          fwdMatches.push_back(new match_s((*C)[2],
                                           S1, (*C)[4], iid1, pos1, len1, ori1,
                                           S2, (*C)[8], iid2, pos2, len2, ori2));
        else
          revMatches.push_back(new match_s((*C)[2],
                                           S1, (*C)[4], iid1, pos1, len1, ori1,
                                           S2, (*C)[8], iid2, pos2, len2, ori2));
      }
    }

    delete C;

    //  Read the next match, if we keep going.  If we are supposed to
    //  stop, we remember inLine for the next call -- it holds the
    //  next match we want!
    //
    if (cont)
      fgets(inLine, 1024, stdin);

  } while (!feof(stdin) && cont);


  fprintf(stderr, " with %8d fwd and %8d rev matches\r", fwdMatches.size(), revMatches.size());
  fflush(stderr);
}





int
main(int argc, char *argv[]) {

  int arg=1;
  while (arg < argc) {

    arg++;
  }

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  char *inLine = new char [1024];
  char *file1  = new char [1024];
  char *file2  = new char [1024];

  readHeader(inLine, stdin, file1, file2, stdout);

  //  cachesize, loadall, report
  //
  FastACache  *C1 = new FastACache(file1, 1, true,  false);
  FastACache  *C2 = new FastACache(file2, 1, false, false);

  vector<match_s *>  fwdMatches;
  vector<match_s *>  revMatches;

  while (!feof(stdin)) {
    readMatches(inLine, C1, C2, fwdMatches, revMatches);

    if (fwdMatches.size() > 0)
      sort(fwdMatches.begin(), fwdMatches.end(), MatchCompare());
    if (revMatches.size() > 0)
      sort(revMatches.begin(), revMatches.end(), MatchCompare());

    u32bit diag_start = 0;
    while (diag_start < fwdMatches.size()) {
#if 0
      fprintf(stderr, "fwd: M u %s . %s %d %d 1 %s %d %d 1\n",
              fwdMatches[diag_start]->_matchId,
              fwdMatches[diag_start]->_id1, fwdMatches[diag_start]->_acc1->getRangeBegin(), fwdMatches[diag_start]->_acc1->getRangeLength(),
              fwdMatches[diag_start]->_id2, fwdMatches[diag_start]->_acc2->getRangeBegin(), fwdMatches[diag_start]->_acc2->getRangeLength());
#endif
      diag_start = extend_matches_on_diagonal(fwdMatches, diag_start);
    }

    diag_start = 0;
    while (diag_start < revMatches.size()) {
#if 0
      fprintf(stderr, "rev: M u %s . %s %d %d 1 %s %d %d 1\n",
              revMatches[diag_start]->_matchId,
              revMatches[diag_start]->_id1, revMatches[diag_start]->_acc1->getRangeBegin(), revMatches[diag_start]->_acc1->getRangeLength(),
              revMatches[diag_start]->_id2, revMatches[diag_start]->_acc2->getRangeBegin(), revMatches[diag_start]->_acc2->getRangeLength());
#endif
      diag_start = extend_matches_on_diagonal(revMatches, diag_start);
    }

    fflush(stdout);
    fflush(stderr);

    //  Dump and destroy all the matches
    //
    for (u32bit i=0; i<fwdMatches.size(); i++) {
      if (!fwdMatches[i]->isDeleted())
        fprintf(stdout, "M u %s . %s "u32bitFMT" "u32bitFMT" 1 %s "u32bitFMT" "u32bitFMT" 1\n",
                fwdMatches[i]->_matchId,
                fwdMatches[i]->_id1, fwdMatches[i]->_acc1->getRangeBegin(), fwdMatches[i]->_acc1->getRangeLength(),
                fwdMatches[i]->_id2, fwdMatches[i]->_acc2->getRangeBegin(), fwdMatches[i]->_acc2->getRangeLength());
      delete fwdMatches[i];
    }

    for (u32bit i=0; i<revMatches.size(); i++) {
      if (!revMatches[i]->isDeleted())
        fprintf(stdout, "M u %s . %s "u32bitFMT" "u32bitFMT" 1 %s "u32bitFMT" "u32bitFMT" -1\n",
                revMatches[i]->_matchId,
                revMatches[i]->_id1, revMatches[i]->_acc1->getRangeBegin(), revMatches[i]->_acc1->getRangeLength(),
                revMatches[i]->_id2, revMatches[i]->_acc2->getRangeBegin(), revMatches[i]->_acc2->getRangeLength());
      delete revMatches[i];
    }

    fflush(stdout);
    fflush(stderr);

    fwdMatches.clear();
    revMatches.clear();
  }
}
