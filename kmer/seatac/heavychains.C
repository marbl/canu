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

  // ignoring id1 and id2 ?
  // Assert that the data is proper:

  if ((xln <= 0) || (yln <= 0) || (xlo < 0) || (ylo < 0)) {
    fprintf(stderr, "StrandPair::addHit()-- bogus data.\n");
    fprintf(stderr, "StrandPair::addHit()-- xln = %d\n", xln);
    fprintf(stderr, "StrandPair::addHit()-- yln = %d\n", yln);
    fprintf(stderr, "StrandPair::addHit()-- xlo = %d\n", xlo);
    fprintf(stderr, "StrandPair::addHit()-- ylo = %d\n", ylo);
    exit(1);
  }

  old_stra1 = id1;
  old_stra2 = id2;

  tmp.ori = direction;
  tmp.xlo = xlo;
  tmp.ylo = ylo;
  tmp.filled = filled;

  // convert to a bounding box
  tmp.xhi = tmp.xlo + xln;
  tmp.yhi = tmp.ylo + yln;
  tmp.xhi = tmp.xlo + xln;
  tmp.yhi = tmp.ylo + yln;
  tmp.S = 0; tmp.neS = 0; tmp.nwS = 0; tmp.seS = 0; tmp.swS = 0;

#if 1
  // Use the match lengths to initialize the self scores.
  if (yln < xln) 
    tmp.selfS = yln;
  else
    tmp.selfS = xln;
#else
  tmp.selfS = filled;
#endif
    
  P.push_back(tmp);
}


  
// new strand pair: begin processing data for the strand pair
//
void StrandPair::process(void) {

  if (beVerbose > 0)
    fprintf(stderr,"HeavyChains: filtering strands %d %d %d\n",old_stra1,old_stra2,P.size());

  if (P.size() > 0) {
    DPTree *dp = NULL;
    dp = new DPTree(P.size(),&(P[0]));
    dp->setParams(maxJump);

    for(int quadrant=0; quadrant < 4; ++quadrant) {
      if (beVerbose > 1)
	fprintf(stderr,"HeavyChains: arranging process quadrant %d\n", quadrant);

      switch(quadrant) {
        case 0:
	case 2:
          for(unsigned i=0; i < P.size(); ++i) {
            int swapi;
            swapi = -P[i].xlo; P[i].xlo = -P[i].xhi; P[i].xhi = swapi;
          }
          break;
        case 1:
        case 3:
          for(unsigned i=0; i < P.size(); ++i) {
            int swapi;
            swapi = -P[i].ylo; P[i].ylo = -P[i].yhi; P[i].yhi = swapi;
          }
          break;
        default:
          fprintf(stderr, "StrandPair::process()-- got invalid quadrant %d\n", quadrant);
          exit(127);
      }

      if (beVerbose > 1)
	fprintf(stderr,"HeavyChains: scoring quadrant\n");
	  
      dp->treeScore();
	  
      if (beVerbose>1) {
	fprintf(stderr,"HeavyChains: recording scores\n");
      }
      switch(quadrant) {
      case 0:
	for(unsigned i=0; i < P.size(); ++i) P[i].nwS = P[i].S;
	break;
      case 1:
	for(unsigned i=0; i < P.size(); ++i) P[i].swS = P[i].S;
	break;
      case 2:
	for(unsigned i=0; i < P.size(); ++i) P[i].seS = P[i].S;
	break;
      case 3:
	for(unsigned i=0; i < P.size(); ++i) P[i].neS = P[i].S;
	break;
      default:
	exit(127);
      }
	  
      if (beVerbose > 1)
	fprintf(stderr,"HeavyChains: done quadrant\n");
    }
    // All output information is now in the match records of P.
    delete dp;
  }
}



long
StrandPair::print(FILE *outF,
                  long matchid) {
  double inc, dec;

  for(unsigned i=0; i < P.size(); ++i) {

    // symmetrize the forward and backward scores
    inc = P[i].neS + P[i].swS - P[i].selfS; // forward complement orientations
    dec = P[i].seS + P[i].nwS - P[i].selfS; // reverse complement orientations

    // Each score already contains the self score 
      
    if ((inc >= minScore) || (dec >= minScore)) {
      int len1 = (P[i].xhi-P[i].xlo);
      int len2 = (P[i].yhi-P[i].ylo);
      matchid++;

      fprintf(outF, "M x H%ld . %s:%d %d %d %d %s:%d %d %d %d > /hf=%.1f /hr=%.1f\n",
              matchid,
              assemblyId1, old_stra1, P[i].xlo, len1, 1,
              assemblyId2, old_stra2, P[i].ylo, len2, (P[i].ori == 'f'? 1 : -1),
              inc, dec );

      sumlen1 += len1;
      sumlen2 += len2;
      maxlen1 = ( maxlen1 > len1 ? maxlen1 : len1);
      maxlen2 = ( maxlen2 > len2 ? maxlen2 : len2);
      maxScoreFwd = (maxScoreFwd > inc ? maxScoreFwd : inc);
      maxScoreRev = (maxScoreRev > dec ? maxScoreRev : dec);
    }
  }

  if (beVerbose > 0)
    fprintf(stderr,"HeavyChains: finished strands %d %d %f %f %f %f\n",
            old_stra1,
            old_stra2,
            maxlen1,
            maxlen2,
            maxScoreFwd,
            maxScoreRev);

  return(matchid);
}
