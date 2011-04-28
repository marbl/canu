//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//


void
cmGlobalData::doSearchRFS(cmComputation *c,
                          cmThreadData  *t) {

  t->pathPos = 0;
  t->pathAdd = 0;

  if (t->path == NULL) {
    t->pathMax = depthMax;
    t->path    = new searchNode [t->pathMax];
  }

  for (uint32 iter=0; iter<500; iter++) {
    t->pathPos = 0;

    t->path[t->pathPos].pIID = c->fragIID;
    t->path[t->pathPos].p5p3 = c->frag5p3;
    t->path[t->pathPos].pLen = fi[c->fragIID].clearLength;
    t->path[t->pathPos].oMax = bbLen[c->fragIID];
    t->path[t->pathPos].oPos = 0;
    t->path[t->pathPos].oLst = bbPos[c->fragIID];

#if 0
    fprintf(stderr, "PATH [%3d] %d/%s' len %d\n",
            t->pathPos,
            t->path[t->pathPos].pIID,
            (t->path[t->pathPos].p5p3 == true) ? "5'3'" : "3'5'",
            t->path[t->pathPos].pLen);
#endif

    //  Follow random paths until we get too long or too deep.  If we find the answer we immediately
    //  return.  If we don't find the answer we exit the while and do another iteration.

    while ((t->path[t->pathPos].pLen < pathMax) &&
           (t->pathPos               < depthMax)) {

      if (testSearch(c, t, tgPos, tgLen))
        //  If any of the target overlaps are the answer
        return;

      //  Compute which edges extend in the correct direction.
      //
      t->extLen = 0;
      for (uint32 test=0; test<bbLen[t->path[t->pathPos].pIID]; test++) {
        overlapInfo  *novl = t->path[t->pathPos].oLst + test;
        uint32        niid = novl->iid;
        bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
        uint32        nlen = 0;

        computeNextPlacement(c, t, novl, niid, n5p3, nlen);

        if (nlen > t->path[t->pathPos].pLen)
          t->ext[t->extLen++] = test;
      }

      if (t->extLen == 0)
        //  If there are no backbone overlaps out of here
        break;

      t->path[t->pathPos].oPos = t->ext[lrand48() % t->extLen];

      overlapInfo  *novl = t->path[t->pathPos].oLst + t->path[t->pathPos].oPos;
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
      uint32        nlen = 0;

      computeNextPlacement(c, t, novl, niid, n5p3, nlen);

      //  We're guaranteed to always advance the path (unlike DFS) so do it.

      t->pathPos++;

      t->path[t->pathPos].pIID = niid;
      t->path[t->pathPos].p5p3 = n5p3;
      t->path[t->pathPos].pLen = nlen;
      t->path[t->pathPos].oMax = bbLen[niid];
      t->path[t->pathPos].oPos = 0;
      t->path[t->pathPos].oLst = bbPos[niid];

      assert(t->path[t->pathPos].pLen > t->path[t->pathPos-1].pLen);

#if 0
      fprintf(stderr, "PATH [%3d] %d/%s' len %d%s\n",
              t->pathPos,
              t->path[t->pathPos].pIID,
              (t->path[t->pathPos].p5p3 == true) ? "5'3'" : "3'5'",
              t->path[t->pathPos].pLen,
              (t->path[t->pathPos].pLen < t->path[t->pathPos-1].pLen) ? "  PATH SHORTER" : "");
#endif
    }

  }  //  Try a bunch of random stabs to find the path

  //  Not found.
  assert(c->sFound == false);
}
