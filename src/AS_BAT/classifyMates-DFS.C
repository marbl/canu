//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//
void
cmGlobalData::doSearchDFS(cmComputation *c,
                          cmThreadData  *t) {
  set<uint32>  visited5p3;
  set<uint32>  visited3p5;

  fprintf(stderr, "SEARCh %d/%s to %d/%s\n",
          c->fragIID, (c->frag5p3 == true) ? "5'3'" : "3'5'", 
          c->mateIID, (c->mate5p3 == true) ? "5'3'" : "3'5'");

  t->pathDepth = 1;  //  CRITICAL; end of loop we increment t->pathPosn[t->pathDepth-1]

  t->pathIID[t->pathDepth]  = c->fragIID;
  t->path5p3[t->pathDepth]  = c->frag5p3;
  t->pathLen[t->pathDepth]  = fi[c->fragIID].clearLength;
  t->pathRoot[t->pathDepth] = bbPos[c->fragIID];
  t->pathPosn[t->pathDepth] = 0;
  t->pathMaxp[t->pathDepth] = bbLen[c->fragIID];

  //  While we still have paths to follow
  //
  while (t->pathDepth > 0) {

    //  Over all overlaps at this depth
    //
    for (;
         t->pathPosn[t->pathDepth] < t->pathMaxp[t->pathDepth];
         t->pathPosn[t->pathDepth]++) {

      //  If we've already seen this fragment in this orientation, get out of here.
      //
      set<uint32> &visited = (t->path5p3[t->pathDepth] == true) ? visited5p3 : visited3p5;

      if (visited.find(t->pathIID[t->pathDepth]) != visited.end())
        continue;

      //  Try extending the fragment at this overlap
      //
      overlapInfo  *novl = t->pathRoot[t->pathDepth] + t->pathPosn[t->pathDepth];
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path5p3[t->pathDepth]) : (t->path5p3[t->pathDepth]);
      uint32        nlen = 0;

      computeNextPlacement(c, t, novl, niid, n5p3, nlen);

      if (nlen < t->pathLen[t->pathDepth])
        //  Went backwards!
        continue;

      //  Went forwards.  Save the extension, unless that would put us over the depth limit.

      t->pathDepth++;

      t->pathIID[t->pathDepth]  = niid;
      t->path5p3[t->pathDepth]  = n5p3;
      t->pathLen[t->pathDepth]  = nlen;
      t->pathRoot[t->pathDepth] = bbPos[niid];
      t->pathPosn[t->pathDepth] = 0;
      t->pathMaxp[t->pathDepth] = bbLen[niid];

      if (testSearch(c, t, tgPos, tgLen))
        //  If any of the target overlaps are the answer
        return;

      if (t->pathDepth ==  depthMax)
        //  End of the line.  Do not pass go.  Proceed directly to, ummm, the next overlap.
        t->pathDepth--;
    }

    //  We've exhausted the paths for this fragment.  Mark it finished and back up one.

    set<uint32> &visited = (t->path5p3[t->pathDepth] == true) ? visited5p3 : visited3p5;
    visited.insert(t->pathIID[t->pathDepth]);

    t->pathDepth--;
  }

  //  Not found.
  assert(c->pathFound == false);
}

