// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Copyright (c) 2005 The J. Craig Venter Institute
// Author: Clark Mobarry
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

#include <errno.h>
#include "heavychains.H"



// The following would need to parameterized for a general kD tree.
// Could we use one function with a static variable to remember the
// sorting direction?
//
int x_compar(const void *x,const void *y) {
  const Match &p1=*((const Match*)x);
  const Match &p2=*((const Match*)y);
  if (p1.xhi < p2.xhi) return -1;
  if (p1.xhi > p2.xhi) return  1;
  return 0;
}

int y_compar(const void *x,const void *y) {
  const Match &p1=*((const Match*)x);
  const Match &p2=*((const Match*)y);
  if (p1.yhi < p2.yhi) return -1;
  if (p1.yhi > p2.yhi) return  1;
  return 0;
}



void StrandPair::addHit(char   direction,
                        u32bit id1,
                        u32bit xlo,
                        u32bit xln,
                        u32bit id2,
                        u32bit ylo,
                        u32bit yln,
                        u32bit filled) {
  Match tmp;

  tmp.xlo = xlo;
  tmp.ylo = ylo;

  tmp.xhi = xlo + xln;
  tmp.yhi = ylo + yln;

  // Use the match lengths to initialize the self scores.
  tmp.selfS = xln;
  if (yln < xln) 
    tmp.selfS = yln;

  tmp.S     = 0.0;
  tmp.neS   = 0;
  tmp.nwS   = 0;
  tmp.seS   = 0;
  tmp.swS   = 0;

  tmp.filled = filled;
  tmp.ori    = direction;

  iid1 = id1;
  iid2 = id2;

  if (beVerbose > 1)
    fprintf(stderr, "heavychains: add %8d %8d %8d -- %8d %8d %8d\n", id1, tmp.xlo, tmp.xhi, id2, tmp.ylo, tmp.yhi);

  Padd(&tmp);
}


  
// new strand pair: begin processing data for the strand pair
//
void StrandPair::process(void) {
  int swapi;

  if (Plen > 0) {
    if (beVerbose > 0)
      fprintf(stderr,"HeavyChains: filtering strands "u32bitFMT" "u32bitFMT" "u32bitFMT"\n", iid1, iid2, Plen);

    DPTree *dp = NULL;
    dp = new DPTree(Plen, P);
    dp->setParams(maxJump);

    for(int quadrant=0; quadrant < 4; ++quadrant) {
      if (beVerbose > 1)
	fprintf(stderr,"HeavyChains: arranging process quadrant %d\n", quadrant);

      if ((quadrant == 0) || (quadrant == 2)) {
        for (int i=0; i<Plen; ++i) {
          swapi    = -P[i].xlo;
          P[i].xlo = -P[i].xhi;
          P[i].xhi = swapi;
        }
      } else {
        for (int i=0; i<Plen; ++i) {
          swapi    = -P[i].ylo;
          P[i].ylo = -P[i].yhi;
          P[i].yhi = swapi;
        }
      }

      if (beVerbose > 1)
	fprintf(stderr,"HeavyChains: scoring quadrant\n");

      dp->treeScore();

      if (beVerbose>1)
	fprintf(stderr,"HeavyChains: recording scores\n");

      switch(quadrant) {
        case 0: for (int i=0; i < Plen; ++i) P[i].nwS = P[i].S; break;
        case 1: for (int i=0; i < Plen; ++i) P[i].swS = P[i].S; break;
        case 2: for (int i=0; i < Plen; ++i) P[i].seS = P[i].S; break;
        case 3: for (int i=0; i < Plen; ++i) P[i].neS = P[i].S; break;
      }
	  
      if (beVerbose > 1)
	fprintf(stderr,"HeavyChains: done quadrant\n");
    }

    // All output information is now in the match records of P.
    delete dp;
  }
}



u64bit
StrandPair::print(FILE   *outF,
                  u64bit  matchid) {

  for (int i=0; i<Plen; ++i) {

    // symmetrize the forward and backward scores
    double inc = P[i].neS + P[i].swS - P[i].selfS; // forward complement orientations
    double dec = P[i].seS + P[i].nwS - P[i].selfS; // reverse complement orientations

    // Each score already contains the self score 
      
    if ((inc >= minScore) || (dec >= minScore)) {
      int len1 = (P[i].xhi-P[i].xlo);
      int len2 = (P[i].yhi-P[i].ylo);
      matchid++;

      if (beVerbose > 1)
        fprintf(stderr, "heavychains: out "u32bitFMTW(8)" %8d %8d -- "u32bitFMTW(8)" %8d %8d\n",
                iid1, P[i].xlo, P[i].xhi,
                iid2, P[i].ylo, P[i].yhi);

      errno = 0;
      fprintf(outF, "M x H"u64bitFMT" . %s:"u32bitFMT" %d %d %d %s:"u32bitFMT" %d %d %d > /hf=%.1f /hr=%.1f\n",
              matchid,
              assemblyId1, iid1, P[i].xlo, len1, 1,
              assemblyId2, iid2, P[i].ylo, len2, (P[i].ori == 'f'? 1 : -1),
              inc, dec);
      if (errno)
        fprintf(stderr, "StrandPair::print()-- write failed: %s\n", strerror(errno));

      sumlen1    += len1;
      sumlen2    += len2;
      maxlen1     = (maxlen1 > len1) ? maxlen1 : len1;
      maxlen2     = (maxlen2 > len2) ? maxlen2 : len2;
      maxScoreFwd = (maxScoreFwd > inc) ? maxScoreFwd : inc;
      maxScoreRev = (maxScoreRev > dec) ? maxScoreRev : dec;
    }

    if (beVerbose > 0)
      fprintf(stderr, "HeavyChains: finished strands "u32bitFMTW(8)" "u32bitFMTW(8)" maxlen1=%f maxlen2=%f maxScoreFwd=%f maxScoreRef=%f\n",
              iid1, iid2, maxlen1, maxlen2, maxScoreFwd, maxScoreRev);
  }

  return(matchid);
}
