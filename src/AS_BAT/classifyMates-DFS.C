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

  t->pathPos = 0;
  t->pathAdd = 0;

  if (t->path == NULL) {
    t->pathMax = depthMax;
    t->path    = new searchNode [t->pathMax];
  }

#if 0
  fprintf(stderr, "SEARCH %d/%s to %d/%s\n",
          c->fragIID, (c->frag5p3 == true) ? "5'3'" : "3'5'", 
          c->mateIID, (c->mate5p3 == true) ? "5'3'" : "3'5'");
#endif

  t->pathPos = 1;  //  CRITICAL; end of loop we increment t->pathOPos[t->pathPos-1]

  t->path[t->pathPos].pIID = c->fragIID;
  t->path[t->pathPos].p5p3 = c->frag5p3;
  t->path[t->pathPos].pLen = fi[c->fragIID].clearLength;
  t->path[t->pathPos].oMax = bbLen[c->fragIID];
  t->path[t->pathPos].oPos = 0;
  t->path[t->pathPos].oLst = bbPos[c->fragIID];

  //  While we still have paths to follow
  //
  while (t->pathPos > 0) {

    //  Over all overlaps at this depth
    //
    for (;
         t->path[t->pathPos].oPos < t->path[t->pathPos].oMax;
         t->path[t->pathPos].oPos++) {

      //  If we've already seen this fragment in this orientation, get out of here.
      //
      set<uint32> &visited = (t->path[t->pathPos].p5p3 == true) ? visited5p3 : visited3p5;

      if (visited.find(t->path[t->pathPos].pIID) != visited.end())
        continue;

      //  Try extending the fragment at this overlap
      //
      overlapInfo  *novl = t->path[t->pathPos].oLst + t->path[t->pathPos].oPos;
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
      uint32        nlen = 0;

      computeNextPlacement(c, t, novl, niid, n5p3, nlen);

      if (nlen < t->path[t->pathPos].pLen)
        //  Went backwards!
        continue;

      //  Went forwards.  Save the extension, unless that would put us over the depth limit.

      t->pathPos++;

      t->path[t->pathPos].pIID  = niid;
      t->path[t->pathPos].p5p3  = n5p3;
      t->path[t->pathPos].pLen  = nlen;
      t->path[t->pathPos].oMax = bbLen[niid];
      t->path[t->pathPos].oPos = 0;
      t->path[t->pathPos].oLst = bbPos[niid];

      if (testSearch(c, t, tgPos, tgLen))
        //  If any of the target overlaps are the answer
        return;

      if (t->pathPos == depthMax)
        //  End of the line.  Do not pass go.  Proceed directly to, ummm, the next overlap.
        t->pathPos--;
    }

    //  We've exhausted the paths for this fragment.  Mark it finished and back up one.

    set<uint32> &visited = (t->path[t->pathPos].p5p3 == true) ? visited5p3 : visited3p5;
    visited.insert(t->path[t->pathPos].pIID);

    t->pathPos--;
  }

  //  Not found.
  assert(c->sFound == false);
}

