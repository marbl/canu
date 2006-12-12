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


u32bit  minEndRunLen = 10;    // -E /matchExtenderMinEndRunLen
u32bit  maxMMBlock   = 3;     // -B /matchExtenderMaxMMBlock
u32bit  minBlockSep  = 20;    // -S /matchExtenderMinBlockSep
double  minIdentity  = 0.95;  // -I /matchExtenderMinIdentity
u32bit  maxNbrSep    = 100;   // -P /matchExtenderMaxNbrSep
u32bit  maxNbrPathMM = 5;     // -D /matchExtenderMaxNbrPathMM


bool   trim_to_pct(vector<match_s *>& matches, u32bit midx, double pct);
void   extend_match_backward(vector<match_s *>& matches, u32bit midx, u32bit min_start_pos);
bool   can_reach_nearby_match(match_s *src, match_s *dest);
bool   extend_match_forward(vector<match_s *>& matches, u32bit midx, match_s *target);
u32bit extend_matches_on_diagonal(vector<match_s *>& matches, u32bit diag_start);


class MatchCompare {
public:
  int operator()(const match_s *m1, const match_s *m2) {
    return(*m1 < *m2);
  }
};



//  Read matches until the iid differs.  Leave the next match in inLine.
//
void
readMatches(atacMatchList     &AL,
            u32bit            &firstMatch,
            FastACache        *C1,
            FastACache        *C2,
            vector<match_s *> &fwdMatches,
            vector<match_s *> &revMatches) {

  fwdMatches.clear();
  revMatches.clear();

  if (firstMatch >= AL.numberOfMatches())
    return;

  u32bit               iid1 = AL.getMatch(firstMatch)->iid1;
  u32bit               iid2 = AL.getMatch(firstMatch)->iid2;

  FastASequenceInCore *seq1 = C1->getSequence(iid1);
  FastASequenceInCore *seq2 = C2->getSequence(iid2);

  while (firstMatch < AL.numberOfMatches()) {
    atacMatch  *m = AL.getMatch(firstMatch++);

    if ((m->iid1 == iid1) && (m->iid2 == iid2)) {
      if (m->fwd1 == m->fwd2)
        fwdMatches.push_back(new match_s(m->matchuid,
                                         seq1, m->iid1, m->pos1, m->len1, m->fwd1,
                                         seq2, m->iid2, m->pos2, m->len2, m->fwd2));
      else
        revMatches.push_back(new match_s(m->matchuid,
                                         seq1, m->iid1, m->pos1, m->len1, m->fwd1,
                                         seq2, m->iid2, m->pos2, m->len2, m->fwd2));
    } else {
      break;
    }
  }

  if (fwdMatches.size() > 0)
    sort(fwdMatches.begin(), fwdMatches.end(), MatchCompare());

  if (revMatches.size() > 0)
    sort(revMatches.begin(), revMatches.end(), MatchCompare());
}




int
main(int argc, char *argv[]) {
  bool    fail = false;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-e") == 0) {
      minEndRunLen = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-b") == 0) {
      maxMMBlock   = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-s") == 0) {
      minBlockSep  = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-i") == 0) {
      minIdentity  = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-p") == 0) {
      maxNbrSep    = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-d") == 0) {
      maxNbrPathMM = strtou32bit(argv[++arg], 0L);
    } else {
      fprintf(stderr, "unknown option %s\n", argv[arg]);
      fail = true;
    }

    arg++;
  }

  if (fail) {
    fprintf(stderr, "usage: %s [options] < matches.atac > matches.atac\n", argv[0]);
    fprintf(stderr, "  -e <int>     matchExtenderMinEndRunLen,  10\n");
    fprintf(stderr, "  -b <int>     matchExtenderMaxMMBlock,     3\n");
    fprintf(stderr, "  -s <int>     matchExtenderMinBlockSep,   20\n");
    fprintf(stderr, "  -i <float>   matchExtenderMinIdentity,    0.95\n");
    fprintf(stderr, "  -p <int>     matchExtenderMaxNbrSep,    100\n");
    fprintf(stderr, "  -d <int>     matchExtenderMaxNbrPathMM,   5\n");
    exit(1);
  }

  atacFile        AF("-");
  atacMatchList  &AM = *AF.matches();

  FastACache  *C1 = new FastACache(AF.assemblyFileA(), 1, true,  false);
  FastACache  *C2 = new FastACache(AF.assemblyFileB(), 1, false, false);

  vector<match_s *>  fwdMatches;
  vector<match_s *>  revMatches;

  u32bit       firstMatch = 0;

  while (firstMatch < AM.numberOfMatches()) {
    readMatches(AM, firstMatch, C1, C2, fwdMatches, revMatches);

    u32bit diag_start = 0;
    while (diag_start < fwdMatches.size()) {
      //fprintf(stderr, "fwd: M u %s . %s %d %d 1 %s %d %d 1\n",
      //        fwdMatches[diag_start]->_matchId,
      //        fwdMatches[diag_start]->_id1, fwdMatches[diag_start]->_acc1->getRangeBegin(), fwdMatches[diag_start]->_acc1->getRangeLength(),
      //        fwdMatches[diag_start]->_id2, fwdMatches[diag_start]->_acc2->getRangeBegin(), fwdMatches[diag_start]->_acc2->getRangeLength());
      diag_start = extend_matches_on_diagonal(fwdMatches, diag_start);
    }

    diag_start = 0;
    while (diag_start < revMatches.size()) {
      //fprintf(stderr, "rev: M u %s . %s %d %d 1 %s %d %d 1\n",
      //        revMatches[diag_start]->_matchId,
      //        revMatches[diag_start]->_id1, revMatches[diag_start]->_acc1->getRangeBegin(), revMatches[diag_start]->_acc1->getRangeLength(),
      //        revMatches[diag_start]->_id2, revMatches[diag_start]->_acc2->getRangeBegin(), revMatches[diag_start]->_acc2->getRangeLength());
      diag_start = extend_matches_on_diagonal(revMatches, diag_start);
    }


    //  Dump and destroy all the matches
    //
    for (u32bit i=0; i<fwdMatches.size(); i++) {
      if (!fwdMatches[i]->isDeleted())
        fprintf(stdout, "M u %s:"u32bitFMT" . %s "u32bitFMT" "u32bitFMT" 1 %s "u32bitFMT" "u32bitFMT" 1\n",
                fwdMatches[i]->_matchId,
                AF.labelA(), fwdMatches[i]->_iid1, fwdMatches[i]->_acc1->getRangeBegin(), fwdMatches[i]->_acc1->getRangeLength(),
                AF.labelB(), fwdMatches[i]->_iid2, fwdMatches[i]->_acc2->getRangeBegin(), fwdMatches[i]->_acc2->getRangeLength());
      delete fwdMatches[i];
    }

    for (u32bit i=0; i<revMatches.size(); i++) {
      if (!revMatches[i]->isDeleted())
        fprintf(stdout, "M u %s:"u32bitFMT" . %s "u32bitFMT" "u32bitFMT" 1 %s "u32bitFMT" "u32bitFMT" -1\n",
                revMatches[i]->_matchId,
                AF.labelA(), revMatches[i]->_iid1, revMatches[i]->_acc1->getRangeBegin(), revMatches[i]->_acc1->getRangeLength(),
                AF.labelB(), revMatches[i]->_iid2, revMatches[i]->_acc2->getRangeBegin(), revMatches[i]->_acc2->getRangeLength());
      delete revMatches[i];
    }
  }
}
