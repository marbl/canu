
class
saveDistance {
public:
  saveDistance() {
  };

  ~saveDistance() {
  };

  //  We save overlaps based on the size of the overlap.
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

#if NDISTD > 0
  int32   doveDist5arr[NDISTD];  //  scratch array for finding the nth largest distance
#endif
#if NDISTC > 0
  int32   coneDist5arr[NDISTC];
  int32   conrDist5arr[NDISTC];
#endif

#if NDISTD > 0
  int32   doveDist3arr[NDISTD];
#endif
#if NDISTC > 0
  int32   coneDist3arr[NDISTC];
  int32   conrDist3arr[NDISTC];
#endif

  int32   doveDist5;     //  minimum distance we should be saving
  int32   coneDist5;
  int32   conrDist5;

  int32   doveDist3;
  int32   coneDist3;
  int32   conrDist3;



  //  Save the N largest values - sorted largest to smallest
#if NDISTD > 0
  void    saveDistMax(int32 *darr, int32  dist) {

    assert(dist >= 0);

    if (dist < darr[NDISTD-1])
      //  Distance less than the smallest distance we want to keep, don't save
      return;

    //  We either want to save a new distance, pushing the last one off of the array,
    //  or notice that we've already saved this distance and leave the array alone.

    for (int32 i=0; i<NDISTD; i++) {
      if (darr[i] == dist)
        //  Saved it already.
        return;

      if (darr[i] < dist) {
        //  Save at i, push i and following down one slot.
        for (int32 j=NDISTD-1; j>i; j--)
          darr[j] = darr[j-1];
        darr[i] = dist;
        return;
      }
    }

    //  Fail, we should never get here.
    assert(0);
  };
#endif



  //  Save the N smallest values - sorted smallest to largest
#if NDISTC > 0
  void    saveDistMin(int32 *darr, int32  dist) {

    assert(dist >= 0);

    if (dist > darr[NDISTC-1])
      //  Distance more than the biggest distance we want to keep, don't save
      return;

    //  We either want to save a new distance, pushing the last one off of the array,
    //  or notice that we've already saved this distance and leave the array alone.

    for (int32 i=0; i<NDISTC; i++) {
      if (darr[i] == dist)
        //  Saved it already.
        return;

      if (darr[i] > dist) {
        //  Save at i, push i and following down one slot.
        for (int32 j=NDISTC-1; j>i; j--)
          darr[j] = darr[j-1];
        darr[i] = dist;
        return;
      }
    }

    //  Fail, we should never get here.
    assert(0);
  };
#endif



  void   compute(fragmentInfo     *fi,
                 OVSoverlap       *ovl,
                 uint32            ovlLen) {

    if (ovlLen == 0)
      return;

#if NDISTD > 0
    for (int32 i=0; i<NDISTD; i++)
      doveDist5arr[i] = doveDist3arr[i] = INT32_MIN;
#endif
#if NDISTC > 0
    for (int32 i=0; i<NDISTC; i++) {
      coneDist5arr[i] = coneDist3arr[i] = INT32_MAX;
      conrDist5arr[i] = conrDist3arr[i] = INT32_MAX;
    }
#endif

    for (uint32 i=0; i<ovlLen; i++) {
      int32  ah = ovl[i].dat.ovl.a_hang;
      int32  bh = ovl[i].dat.ovl.b_hang;
      int32  fa = fi[ovl[i].a_iid].clearLength;
      int32  fb = fi[ovl[i].b_iid].clearLength;

#if NDISTD > 0
      if        (AS_OVS_overlapAEndIs5prime(ovl[i])) {
        //  ah < 0 && bh < 0
        saveDistMax(doveDist5arr, fb - -ah);

      } else if (AS_OVS_overlapAEndIs3prime(ovl[i])) {
        //  ah > 0 && bh > 0
        saveDistMax(doveDist3arr, fb -  bh);
      }
#endif

#if NDISTC > 0
      } else if (AS_OVS_overlapAIsContained(ovl[i])) {
        //  ah <= 0 && bh >= 0
        saveDistMin(coneDist5arr, -ah);
        saveDistMin(coneDist3arr,  bh);
        
      } else if (AS_OVS_overlapAIsContainer(ovl[i])) {
        //  ah >= 0 && bh <= 0
        saveDistMin(conrDist5arr,  ah);
        saveDistMin(conrDist3arr, -bh);

      } else {
        assert(0);
      }
#endif
    }

#if NDISTD > 0
    doveDist5 = doveDist5arr[NDISTD-1];
#endif
#if NDISTC > 0
    coneDist5 = coneDist5arr[NDISTC-1];
    conrDist5 = conrDist5arr[NDISTC-1];
#endif

#if NDISTD > 0
    doveDist3 = doveDist3arr[NDISTD-1];
#endif
#if NDISTC > 0
    coneDist3 = coneDist3arr[NDISTC-1];
    conrDist3 = conrDist3arr[NDISTC-1];
#endif
  }
};



