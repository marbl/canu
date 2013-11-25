
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. Craig Venter Institute
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

#ifndef AS_ALN_BRUTEFORCEDP
#define AS_ALN_BRUTEFORCEDP

static const char *rcsid_AS_ALN_BRUTEFORCEDP = "$Id$";


typedef struct {
  unsigned int  score  : 30;
  unsigned int  action : 2;
} dpCell;



class dpActions {
public:
  dpActions() {
    actions = new uint64 * [2 * AS_READ_MAX_NORMAL_LEN];
    memset(actions, 0, 2 * AS_READ_MAX_NORMAL_LEN * sizeof(uint64 *));
  };

  ~dpActions() {
    for (uint32 i=0; i<2 * AS_READ_MAX_NORMAL_LEN; i++)
      if (actions[i] != NULL)
        delete [] actions[i];

    delete [] actions;
  };

  //  Pointer to array of packed bits.
  uint64  **actions;
  uint32   *actionLength;

  void      allocate(uint32 i) {
    if (actions[i] != NULL)
      return;

    uint32  len = 2 * AS_READ_MAX_NORMAL_LEN * 2 / 64 + 1;

    actions[i] = new uint64 [len];
    memset(actions[i], 0, len * sizeof(uint64));  //  IMPORTANT, init to STOP
  };

  uint64    get(uint32 i, uint32 j) {
    allocate(i);

    uint32   addr = ((j >> 5) & 0x07ffffff);
    uint32   bits = ((j     ) & 0x0000001f) << 1;

    uint64   v = (actions[i][addr] >> bits) & 0x00000003llu;

    //fprintf(stderr, "GET %u %u  %lu (%u %u)\n", i, j, v, addr, bits);

    return(v);
  };

  void      set(uint32 i, uint32 j, uint64 v) {
    allocate(i);

    assert(v < 4);

    uint32   addr = ((j >> 5) & 0x07ffffff);
    uint32   bits = ((j     ) & 0x0000001f) << 1;

    //fprintf(stderr, "SET %u %u to %lu (%u %u)\n", i, j, v, addr, bits);

    actions[i][addr] &= ~(0x00000003llu << bits);
    actions[i][addr] |=  (v             << bits);
  };
};



typedef struct {
  int32   matches;
  int32   alignLen;
  int32   begI, begJ;
  int32   endI, endJ;
  int32   lenA, lenB;

  double  pIdentity;
  double  pCoverageA;
  double  pCoverageB;
} alignLinker_s;


void
alignLinker(char           *alignA,
            char           *alignB,
            char           *stringA,
            char           *stringB,
            dpCell        (*M)[AS_READ_MAX_NORMAL_LEN + 1],
            alignLinker_s  *a,
            int             endToEnd,
            int             allowNs,
            int             ahang, int bhang);

#endif  //  AS_ALN_BRUTEFORCEDP
