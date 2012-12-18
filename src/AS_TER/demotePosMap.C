
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: demotePosMap.C,v 1.1 2012-12-18 07:25:01 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "AS_UTL_splitToWords.H"

#include <string>
#include <vector>
#include <set>
#include <map>

using namespace std;

//  Given a frgscf and ctgscf, build frgctg.
//
//  Expected   promotePosMap ctgscf < frgscf > frgctg
//
//
//  -----
//    -------
//               --------
//                   --------- FRG
//  ---------....------------- SCF
//  ---------    ------------- CTG

class contigDat {
public:
  contigDat(char *name_, uint32 bgn_, uint32 end_, char ori_) {
    memset(name, 0, sizeof(char) * 256);
    strcpy(name, name_);
    bgn = bgn_;
    end = end_;
    ori = ori_;
  };
  ~contigDat() {
  };

  char    name[256];

  uint32  bgn;
  uint32  end;
  char    ori;
};


class scaffoldDat {
public:
  scaffoldDat(char *name_) {
    memset(name, 0, sizeof(char) * 256);
    strcpy(name, name_);

    scfToCtgLen = 8192;
    scfToCtg    = new uint32 [scfToCtgLen];

    for (uint32 i=0; i<scfToCtgLen; i++)
      scfToCtg[i] = UINT32_MAX;
  };
  ~scaffoldDat() {
  };

  void   add(uint32 ctgId, uint32 bgn, uint32 end) {
    while (scfToCtgLen <= end) {
      //fprintf(stderr, "REALLOCATE from %u to %u\n", scfToCtgLen, scfToCtgLen*2);

      uint32 *N = new uint32 [scfToCtgLen * 2];
      uint32  i;

      for (i=0; i<scfToCtgLen; i++)
        N[i] = scfToCtg[i];

      scfToCtgLen *= 2;

      for (; i<scfToCtgLen; i++)
        N[i] = UINT32_MAX;

      delete [] scfToCtg;
      scfToCtg = N;
    }

    //fprintf(stderr, "ADD() ctg %u from %u to %u\n", ctgId, bgn, end);

    for (uint32 i=bgn; i<=end; i++) {
      if (scfToCtg[i] != UINT32_MAX)
        fprintf(stderr, "i=%u scfToCtg[i]=%u\n", i, scfToCtg[i]);
      assert(scfToCtg[i] == UINT32_MAX);
      scfToCtg[i] = ctgId;
    }
  };

  uint32    translate(uint32 coord) {
    if (coord >= scfToCtgLen)
      coord = scfToCtgLen - 1;

    if (scfToCtg[coord] != UINT32_MAX)
      return(scfToCtg[coord]);

    //  No valid value.  Pick the closest.

    int32   negVal=UINT32_MAX, neg=coord-1;
    int32   posVal=UINT32_MAX, pos=coord+1;

    while ((neg > 0) && (scfToCtg[neg] == UINT32_MAX))
      neg--;

    while ((pos < scfToCtgLen-1) && (scfToCtg[pos] == UINT32_MAX))
      pos++;

    fprintf(stderr, "NO VALID for coord %u -- neg %u %u pos %u %u\n",
            coord, neg, scfToCtg[neg], pos, scfToCtg[pos]);

    if (scfToCtg[neg] == UINT32_MAX)
      return(scfToCtg[pos]);

    if (scfToCtg[pos] == UINT32_MAX)
      return(scfToCtg[neg]);

    if (coord-neg > pos-coord)
      return(scfToCtg[pos]);
    else
      return(scfToCtg[neg]);
  };


  char            name[256];

  uint32          bgn;
  uint32          end;

  uint32          scfToCtgLen;
  uint32         *scfToCtg;
};



map<string,uint32>   nameToCtgId;
map<string,uint32>   nameToScfId;

vector<contigDat>    ctgDat;
vector<scaffoldDat>  scfDat;


int
main(int argc, char **argv) {

  //  Load contig to scaffold map

  errno = 0;
  FILE *ctgscf = fopen(argv[1], "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", argv[1], strerror(errno)), exit(1);

  char          L[1024 * 1024];
  splitToWords  W;

  memset(L, 0, sizeof(char) * 1024 * 1024);

  fgets(L, 1024 * 1024, ctgscf);
  assert(L[1024 * 1024 - 1] == 0);

  while (!feof(ctgscf)) {
    W.split(L);

    string  ctgN(W[0]);
    string  scfN(W[1]);

    uint32  ctgId;
    uint32  scfId;

    uint32  bgn = W(2);
    uint32  end = W(3);

    char    ori = W[4][0];

    if (nameToCtgId.count(ctgN) == 0) {
      ctgId = nameToCtgId[ctgN] = ctgDat.size();
      ctgDat.push_back(contigDat(W[0], bgn, end, ori));
    } else {
      ctgId = nameToCtgId[ctgN];
    }

    if (nameToScfId.count(scfN) == 0) {
      scfId = nameToScfId[scfN] = scfDat.size();
      scfDat.push_back(scaffoldDat(W[1]));
    } else {
      scfId = nameToScfId[scfN];
    }

    contigDat    *ctg = &ctgDat[ctgId];
    scaffoldDat  *scf = &scfDat[scfId];

    scf->add(ctgId, bgn, end);

    fgets(L, 1024 * 1024, ctgscf);
    assert(L[1024 * 1024 - 1] == 0);
  }

  fclose(ctgscf);


  //  Do the promotion

  fgets(L, 1024 * 1024, stdin);
  assert(L[1024 * 1024 - 1] == 0);

  while (!feof(stdin)) {
    W.split(L);

    string  frgN(W[0]);
    string  scfN(W[1]);

    uint32  scfId;

    uint32  scfbgn = W(2);
    uint32  scfend = W(3);

    char    scfori = W[4][0];

    if (nameToScfId.count(scfN) == 0) {
      scfId = nameToScfId[scfN] = scfDat.size();
      scfDat.push_back(scaffoldDat(W[1]));
    } else {
      scfId = nameToScfId[scfN];
    }

    scaffoldDat  *scf = &scfDat[scfId];

    //  CHECK FOR OUT OF BOUNDS
    uint32  ctgId1 = scf->translate(scfbgn);
    uint32  ctgId2 = scf->translate(scfend);

    if (ctgId1 != ctgId2)
      fprintf(stderr, "CTG MISMATCH for %s %s %u %u %c\n", W[0], W[1], scfbgn, scfend, scfori);

    contigDat *ctg = &ctgDat[ctgId1];

    uint32  ctgbgn = scfbgn - ctg->bgn;
    uint32  ctgend = scfend - ctg->bgn;
    uint32  ctgori = (scfori == 'f') ? 'f' : 'r';

    if (ctg->ori == 'r') {
      ctgbgn = ctg->end - scfend;
      ctgend = ctg->end - scfbgn;
      ctgori = (scfori == 'f') ? 'r' : 'f';
    }

    fprintf(stdout, "%s\t%s\t%u\t%u\t%c\n",
            W[0], ctg->name, ctgbgn, ctgend, ctgori);

    fgets(L, 1024 * 1024, stdin);
    assert(L[1024 * 1024 - 1] == 0);
  }
}
