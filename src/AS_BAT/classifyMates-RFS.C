//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//



void
cmGlobalData::doSearchRFS(cmComputation *c,
                          cmThreadData  *t) {
  bool         pathFound = false;
  set<uint32>  visited5p3;
  set<uint32>  visited3p5;

  //  If pathInnie == false, then we're attempting to find a path for outtie oriented fragments.
  //  In this case, we start with the first fragment 5p3=false, and need to end with 5p3=true.

  bool    bgn5p3 = (pathInnie == false) ? false : true;
  bool    end5p3 = (pathInnie == false) ? true : false;

  for (uint32 iter=0; iter<850; iter++) {
    c->pathDepth                  = 1;

    c->pathIID[c->pathDepth]  = c->iid;
    c->path5p3[c->pathDepth]  = bgn5p3;
    c->pathLen[c->pathDepth]  = fi[c->iid].clearLength;
    c->pathRoot[c->pathDepth] = oiPos[c->iid];
    c->pathPosn[c->pathDepth] = 0;
    c->pathMaxp[c->pathDepth] = oiLen[c->iid];

    //  Follow random paths until we get too long or too deep.  If we find the answer we immediately
    //  return.  If we don't find the answer we exit the while and do another iteration.

    while ((c->pathLen[c->pathDepth] < pathMax) &&
           (c->pathDepth             < depthMax)) {

      //  Test if any of our overlaps are the answer we're looking for

      for (uint32 test=0; test<c->pathMaxp[c->pathDepth]; test++) {
        overlapInfo  *novl = c->pathRoot[c->pathDepth] + test;

        uint32        niid = novl->iid;                                                         //  Next fragment
        bool          n5p3 = (novl->flipped) ? (!c->path5p3[c->pathDepth]) : (c->path5p3[c->pathDepth]);    //  Next fragment orientation;
        uint32        nlen = 0;                                                                 //  New length of the path, if zero, can't extend

        computeNextPlacement(c, novl, niid, n5p3, nlen);

        if ((niid == fi[c->iid].mateIID) &&
            (n5p3 == end5p3) &&
            (nlen >= pathMin) &&
            (nlen <= pathMax)) {
          c->pathDepth++;

          c->pathIID[c->pathDepth]  = niid;
          c->path5p3[c->pathDepth]  = n5p3;
          c->pathLen[c->pathDepth]  = nlen;
          c->pathRoot[c->pathDepth] = oiPos[niid];
          c->pathPosn[c->pathDepth] = 0;
          c->pathMaxp[c->pathDepth] = oiLen[niid];

          c->pathFound = true;
          return;
        }
      }
      
      //  Pick a random edge out of here and follow it

      if (c->pathMaxp[c->pathDepth] == 0)
        break;

      c->pathPosn[c->pathDepth] = lrand48() % c->pathMaxp[c->pathDepth];

      overlapInfo  *novl = c->pathRoot[c->pathDepth] + c->pathPosn[c->pathDepth];

      uint32        niid = novl->iid;                                                         //  Next fragment
      bool          n5p3 = (novl->flipped) ? (!c->path5p3[c->pathDepth]) : (c->path5p3[c->pathDepth]);    //  Next fragment orientation;
      uint32        nlen = 0;                                                                 //  New length of the path, if zero, can't extend

      computeNextPlacement(c, novl, niid, n5p3, nlen);

      //
      //  If this is an extension, go into it, otherwise stop the search.
      //

      //if (nlen > c->pathLen[c->pathDepth]) {
        c->pathDepth++;

        c->pathIID[c->pathDepth]  = niid;
        c->path5p3[c->pathDepth]  = n5p3;
        c->pathLen[c->pathDepth]  = nlen;
        c->pathRoot[c->pathDepth] = oiPos[niid];
        c->pathPosn[c->pathDepth] = 0;
        c->pathMaxp[c->pathDepth] = oiLen[niid];
        //} else {
        //  break;
        //}
    }

  }  //  Try a bunch of random stabs to find the path

  c->pathFound = false;
  return;
}
