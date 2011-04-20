//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//


void
cmGlobalData::doSearchRFS(cmComputation *c,
                          cmThreadData  *t) {

  //  c is fully initialized by cmComputation()

  for (uint32 iter=0; iter<500; iter++) {
    t->pathDepth = 0;

    t->pathIID[t->pathDepth]  = c->fragIID;
    t->path5p3[t->pathDepth]  = c->frag5p3;
    t->pathLen[t->pathDepth]  = fi[c->fragIID].clearLength;
    t->pathRoot[t->pathDepth] = bbPos[c->fragIID];
    t->pathPosn[t->pathDepth] = 0;
    t->pathMaxp[t->pathDepth] = bbLen[c->fragIID];

#if 0
    fprintf(stderr, "PATH [%3d] %d/%s' len %d\n",
            t->pathDepth,
            t->pathIID[t->pathDepth],
            (t->path5p3[t->pathDepth] == true) ? "5'3'" : "3'5'",
            t->pathLen[t->pathDepth]);
#endif

    //  Follow random paths until we get too long or too deep.  If we find the answer we immediately
    //  return.  If we don't find the answer we exit the while and do another iteration.

    while ((t->pathLen[t->pathDepth] < pathMax) &&
           (t->pathDepth             < depthMax)) {

      if (testSearch(c, t, tgPos, tgLen))
        //  If any of the target overlaps are the answer
        return;

      //  Compute which edges extend in the correct direction.
      //
      t->extLen = 0;
      for (uint32 test=0; test<bbLen[t->pathIID[t->pathDepth]]; test++) {
        overlapInfo  *novl = t->pathRoot[t->pathDepth] + test;
        uint32        niid = novl->iid;
        bool          n5p3 = (novl->flipped) ? (!t->path5p3[t->pathDepth]) : (t->path5p3[t->pathDepth]);
        uint32        nlen = 0;

        computeNextPlacement(c, t, novl, niid, n5p3, nlen);

        if (nlen > t->pathLen[t->pathDepth])
          t->ext[t->extLen++] = test;
      }

      if (t->extLen == 0)
        //  If there are no backbone overlaps out of here
        break;

      t->pathPosn[t->pathDepth] = t->ext[lrand48() % t->extLen];

      overlapInfo  *novl = t->pathRoot[t->pathDepth] + t->pathPosn[t->pathDepth];
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path5p3[t->pathDepth]) : (t->path5p3[t->pathDepth]);
      uint32        nlen = 0;

      computeNextPlacement(c, t, novl, niid, n5p3, nlen);

      //  We're guaranteed to always advance the path (unlike DFS) so do it.

      t->pathDepth++;

      t->pathIID[t->pathDepth]  = niid;
      t->path5p3[t->pathDepth]  = n5p3;
      t->pathLen[t->pathDepth]  = nlen;
      t->pathRoot[t->pathDepth] = bbPos[niid];
      t->pathPosn[t->pathDepth] = 0;
      t->pathMaxp[t->pathDepth] = bbLen[niid];

      assert(t->pathLen[t->pathDepth] > t->pathLen[t->pathDepth-1]);

#if 0
      fprintf(stderr, "PATH [%3d] %d/%s' len %d%s\n",
              t->pathDepth,
              t->pathIID[t->pathDepth],
              (t->path5p3[t->pathDepth] == true) ? "5'3'" : "3'5'",
              t->pathLen[t->pathDepth],
              (t->pathLen[t->pathDepth] < t->pathLen[t->pathDepth-1]) ? "  PATH SHORTER" : "");
#endif
    }

  }  //  Try a bunch of random stabs to find the path

  //  Not found.
  assert(c->pathFound == false);
}
