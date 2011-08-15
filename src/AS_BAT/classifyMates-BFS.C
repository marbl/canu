//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//
void
cmGlobalData::doSearchBFS(cmComputation *c,
                          cmThreadData  *t) {
  set<uint32>               visited5p3;
  set<uint32>               visited3p5;
  map<uint32,overlapInfo*>  solution;

  //  Allocate nodes for our search.  This lets us search 32 million fragments, and takes 1GB of
  //  memory.  We don't try to do anything fancy like round-robin.  Once we full up the list,
  //  we're done.

  t->pathPos    = 0;
  t->pathAdd    = 0;
  t->searchIter = 0;

  if (t->path == NULL) {
    t->pathMax = nodesMax;
    t->path    = new searchNode [t->pathMax];
  }

  //  Build the map from backbone fragment to solution overlap.
  //
  //  The map goes from a-fragID to overlap with the mateIID b-frag.  NOTE that we
  //  DO NOT reset the 'iid' of this overlap to be correct.
  //
  {
    overlapInfo  *pos = gtPos[c->mateIID];
    uint32        len = gtLen[c->mateIID];

    for (uint32 ii=0; ii<len; ii++)
      solution[ pos[ii].iid ] = pos + ii;
  }

  //  Seed the search with the first fragment

  t->path[t->pathAdd].pIID = c->fragIID;
  t->path[t->pathAdd].p5p3 = c->frag5p3;
  t->path[t->pathAdd].pLen = fi[c->fragIID].clearLength;
  t->path[t->pathAdd].oMax = 0;
  t->path[t->pathAdd].oPos = 0;
  t->path[t->pathAdd].oLst = bbPos[c->fragIID];

  t->pathAdd++;

  for (;
       ((t->pathPos < t->pathMax) &&
        (t->pathPos < t->pathAdd));
       t->pathPos++) {

    if ((distMin                  <= t->path[t->pathPos].pLen) &&
        (t->path[t->pathPos].pLen <= distMax) &&
        (testSearch(c, t, solution)))  //  tgPos, tgLen for the old slow method
      //  If any of the target overlaps are the answer
      return;

    if (t->pathAdd >= t->pathMax)
      //  No space for more overlaps, abort adding any.
      continue;

    //  Add more fragments to the search.

    for (uint32 o=0; o < bbLen[t->path[t->pathPos].pIID]; o++) {
      overlapInfo  *novl = t->path[t->pathPos].oLst + o;
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
      uint32        nlen = 0;

      set<uint32> &visited = (t->path[t->pathPos].p5p3 == true) ? visited5p3 : visited3p5;

      if (visited.find(niid) != visited.end())
        //  Been here already.
        continue;

      if (fi[niid].isBackbone == false)
        //  Not a backbone read
        continue;

      visited.insert(niid);

      computeNextPlacement(c, t, novl, niid, n5p3, nlen);

      if (nlen <= t->path[t->pathPos].pLen)
        //  Path went backwards.
        continue;

      if (nlen > distMax)
        //  Path too far, don't add
        continue;

      if (t->pathAdd >= t->pathMax)
        //  No space, don't add
        continue;

      t->path[t->pathAdd].pIID = niid;
      t->path[t->pathAdd].p5p3 = n5p3;
      t->path[t->pathAdd].pLen = nlen;
      t->path[t->pathAdd].oMax = 0;
      t->path[t->pathAdd].oPos = 0;
      t->path[t->pathAdd].oLst = bbPos[niid];

      t->pathAdd++;
      t->searchIter++;
    }
  }

  if (t->pathAdd >= t->pathMax)
    c->sLimited = true;
  else
    c->sExhausted = true;

  //  Not found.
  assert(c->sFound == false);
}

