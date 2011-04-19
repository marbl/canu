//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//


void
cmGlobalData::doSearchRFS(cmComputation *c,
                          cmThreadData  *t) {

  //  If pathInnie == false, then we're attempting to find a path for outtie oriented fragments.
  //  In this case, we start with the first fragment 5p3=false, and need to end with 5p3=true.

  bool    bgn5p3 = (pathInnie == false) ? false : true;
  bool    end5p3 = (pathInnie == false) ? true : false;

  for (uint32 iter=0; iter<500; iter++) {
    c->pathDepth = 0;

    c->pathIID[c->pathDepth]  = c->iid;
    c->path5p3[c->pathDepth]  = bgn5p3;
    c->pathLen[c->pathDepth]  = fi[c->iid].clearLength;
    c->pathRoot[c->pathDepth] = bbPos[c->iid];
    c->pathPosn[c->pathDepth] = 0;
    c->pathMaxp[c->pathDepth] = bbLen[c->iid];

#if 0
    fprintf(stderr, "PATH [%3d] %d/%s' len %d\n",
            c->pathDepth,
            c->pathIID[c->pathDepth],
            (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'",
            c->pathLen[c->pathDepth]);
#endif

    //  Follow random paths until we get too long or too deep.  If we find the answer we immediately
    //  return.  If we don't find the answer we exit the while and do another iteration.

    while ((c->pathLen[c->pathDepth] < pathMax) &&
           (c->pathDepth             < depthMax)) {

      if (testSearch(c, tgPos, tgLen, end5p3))
        //  If any of the target overlaps are the answer
        return;

      //  Compute which edges extend in the correct direction.

      c->extLen = 0;
      for (uint32 test=0; test<bbLen[c->pathIID[c->pathDepth]]; test++) {
        overlapInfo  *novl = c->pathRoot[c->pathDepth] + test;
        uint32        niid = novl->iid;                                                         //  Next fragment
        bool          n5p3 = (novl->flipped) ? (!c->path5p3[c->pathDepth]) : (c->path5p3[c->pathDepth]);    //  Next fragment orientation;
        uint32        nlen = 0;                                                                 //  New length of the path, if zero, can't extend

        computeNextPlacement(c, novl, niid, n5p3, nlen);

        if (nlen > c->pathLen[c->pathDepth])
          c->ext[c->extLen++] = test;
      }

      if (c->extLen == 0)
        //  If there are no backbone overlaps out of here
        break;

      c->pathPosn[c->pathDepth] = c->ext[lrand48() % c->extLen];

      overlapInfo  *novl = c->pathRoot[c->pathDepth] + c->pathPosn[c->pathDepth];
      uint32        niid = novl->iid;                                                         //  Next fragment
      bool          n5p3 = (novl->flipped) ? (!c->path5p3[c->pathDepth]) : (c->path5p3[c->pathDepth]);    //  Next fragment orientation;
      uint32        nlen = 0;                                                                 //  New length of the path, if zero, can't extend

      computeNextPlacement(c, novl, niid, n5p3, nlen);

      //
      //  If this is an extension, go into it, otherwise stop the search.
      //

      c->pathDepth++;

      c->pathIID[c->pathDepth]  = niid;
      c->path5p3[c->pathDepth]  = n5p3;
      c->pathLen[c->pathDepth]  = nlen;
      c->pathRoot[c->pathDepth] = bbPos[niid];
      c->pathPosn[c->pathDepth] = 0;
      c->pathMaxp[c->pathDepth] = bbLen[niid];

      assert(c->pathLen[c->pathDepth] > c->pathLen[c->pathDepth-1]);

#if 0
      fprintf(stderr, "PATH [%3d] %d/%s' len %d%s\n",
              c->pathDepth,
              c->pathIID[c->pathDepth],
              (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'",
              c->pathLen[c->pathDepth],
              (c->pathLen[c->pathDepth] < c->pathLen[c->pathDepth-1]) ? "  PATH SHORTER" : "");
#endif
    }

  }  //  Try a bunch of random stabs to find the path

  c->pathFound = false;
  return;
}
