
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

static const char *rcsid_CLASSIFYMATES_SAVEDISTANCE_H = "$Id: classifyMates-saveDistance.H,v 1.2 2012-02-13 11:41:03 brianwalenz Exp $";

//  Save the longest N overlaps.  Decreasing this value makes the search much much faster, but costs
//  a bit in sensitivity.
//
//  Keep the best and near best overlaps for dovetail overlaps.
//
//    ------------------------
//       ----------------------------- (best)
//          ---------------------- (near best)
//                 --------------
//                     ----------------------------------
//                          -----------------------------------
//
//  Keep containment overlaps that are near the end.
//
//                     ---------
//    ----------------------------- (near the end)
//        -------------------------------
//            ----------------------------------
//                 ------------------------------------ (near the end)
//                     --------------------------------------- (near the end)
//
//  doveDist - Dovetail overlaps with overlap length at least this big are saved.
//             The 'near best' above might not be informative, but it is still
//             the second best overlap and is kept.
//
//  coneDist - Containee overlaps (A contained in B) are saved if they are
//             at least this close to the end of the other fragment.
//
//  conrDist - Container overlaps (A contains B), likewise.

class
saveDistance {
public:
  saveDistance(int32 ndist_) {
    ndist = ndist_;

    doveDist5 = 0;
    doveDist3 = 0;

    doveDist5arr = new int32 [ndist];
    doveDist3arr = new int32 [ndist];
  };

  ~saveDistance() {
    delete [] doveDist5arr;
    delete [] doveDist3arr;
  };

  int32   ndist;

  int32   doveDist5;     //  minimum distance we should be saving
  int32   doveDist3;

  int32  *doveDist5arr;  //  scratch array for finding the nth largest distance
  int32  *doveDist3arr;


  //  Save the N largest values - sorted largest to smallest
  void    saveDistMax(int32 *darr, int32  dist) {

    assert(dist >= 0);

    if (dist < darr[ndist-1])
      //  Distance less than the smallest distance we want to keep, don't save
      return;

    //  We either want to save a new distance, pushing the last one off of the array,
    //  or notice that we've already saved this distance and leave the array alone.

    for (int32 i=0; i<ndist; i++) {
      if (darr[i] == dist)
        //  Saved it already.
        return;

      if (darr[i] < dist) {
        //  Save at i, push i and following down one slot.
        for (int32 j=ndist-1; j>i; j--)
          darr[j] = darr[j-1];
        darr[i] = dist;
        return;
      }
    }

    //  Fail, we should never get here.
    assert(0);
  };


  void   compute(fragmentInfo     *fi,
                 OVSoverlap       *ovl,
                 uint32            ovlLen) {

    if ((ovlLen == 0) || (ndist == 0))
      return;

    memset(doveDist5arr, 0, sizeof(int32) * ndist);
    memset(doveDist3arr, 0, sizeof(int32) * ndist);

    for (uint32 i=0; i<ovlLen; i++) {
      int32  ah = ovl[i].dat.ovl.a_hang;
      int32  bh = ovl[i].dat.ovl.b_hang;
      int32  fa = fi[ovl[i].a_iid].clearLength;
      int32  fb = fi[ovl[i].b_iid].clearLength;

      if        (AS_OVS_overlapAEndIs5prime(ovl[i])) {
        //  ah < 0 && bh < 0
        saveDistMax(doveDist5arr, fb - -ah);

      } else if (AS_OVS_overlapAEndIs3prime(ovl[i])) {
        //  ah > 0 && bh > 0
        saveDistMax(doveDist3arr, fb -  bh);
      }
    }

    doveDist5 = doveDist5arr[ndist-1];
    doveDist3 = doveDist3arr[ndist-1];
  }
};

