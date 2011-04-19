//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//
void
cmGlobalData::doSearchDFS(cmComputation *c,
                          cmThreadData  *t) {
  bool         pathFound = false;
  set<uint32>  visited5p3;
  set<uint32>  visited3p5;

  //  If pathInnie == false, then we're attempting to find a path for outtie oriented fragments.
  //  In this case, we start with the first fragment 5p3=false, and need to end with 5p3=true.

  bool    bgn5p3 = (pathInnie == false) ? false : true;
  bool    end5p3 = (pathInnie == false) ? true : false;

  c->pathDepth                  = 1;

  c->pathIID[c->pathDepth]  = c->iid;
  c->path5p3[c->pathDepth]  = bgn5p3;
  c->pathLen[c->pathDepth]  = fi[c->iid].clearLength;
  c->pathRoot[c->pathDepth] = bbPos[c->iid];
  c->pathPosn[c->pathDepth] = 0;
  c->pathMaxp[c->pathDepth] = bbLen[c->iid];

  //  While we still have paths to follow
  //
  while (c->pathDepth > 0) {

    //  Over all overlaps at this depth
    //
    while (c->pathPosn[c->pathDepth] < c->pathMaxp[c->pathDepth]) {

      //
      //  If we've already seen this fragment in this orientation, get out of here.
      //

      set<uint32> &visited = (c->path5p3[c->pathDepth] == true) ? visited5p3 : visited3p5;

      if (visited.find(c->pathIID[c->pathDepth]) != visited.end()) {
        c->pathPosn[c->pathDepth]++;
        continue;
      }

      //
      //  Are we finished??
      //

      if ((c->pathIID[c->pathDepth] == fi[c->iid].mateIID) &&
          (c->path5p3[c->pathDepth] == end5p3) &&
          (c->pathLen[c->pathDepth] >= pathMin) &&
          (c->pathLen[c->pathDepth] <= pathMax)) {
        c->pathFound = true;
        return;
      }

      //
      //  Try extending into the fragment at this overlap
      //

      overlapInfo  *novl = c->pathRoot[c->pathDepth] + c->pathPosn[c->pathDepth];

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
        c->pathRoot[c->pathDepth] = bbPos[niid];
        c->pathPosn[c->pathDepth] = 0;
        c->pathMaxp[c->pathDepth] = bbLen[niid];

        c->pathFound = true;
        return;
      }

      //
      //  If the length is not too large, extend into it, otherwise, try the next fragment.
      //

      if ((nlen      >  c->pathLen[c->pathDepth]) &&
          (nlen      <= pathMax) &&
          (c->pathDepth <= depthMax)) {
        c->pathDepth++;

        c->pathIID[c->pathDepth]  = niid;
        c->path5p3[c->pathDepth]  = n5p3;
        c->pathLen[c->pathDepth]  = nlen;
        c->pathRoot[c->pathDepth] = bbPos[niid];
        c->pathPosn[c->pathDepth] = 0;
        c->pathMaxp[c->pathDepth] = bbLen[niid];

        //if (VERBOSE4)
        //  fprintf(stderr, "WALK %5u/%s with %5u overlaps at depth %2u of length %5d.\n",
        //          c->pathIID[c->pathDepth], (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'",
        //          c->pathMaxp[c->pathDepth], c->pathDepth, c->pathLen[c->pathDepth]);

        continue;
      }

      //if (VERBOSE5)
      //  fprintf(stderr, "ABRT %5u/%s\n", c->pathIID[c->pathDepth], (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'");

      //  Move to the next overlap for this fragment.
      c->pathPosn[c->pathDepth]++;
    }

    //  We've exhausted the paths for this fragment.  Mark it finished and back up one.

    set<uint32> &visited = (c->path5p3[c->pathDepth] == true) ? visited5p3 : visited3p5;
    visited.insert(c->pathIID[c->pathDepth]);

    c->pathDepth--;
    c->pathPosn[c->pathDepth]++;
  }

  c->pathFound = false;
  return;
}

