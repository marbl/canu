//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//
void
cmGlobalData::doSearchBFS(cmComputation *c,
                          cmThreadData  *t) {

  //  NOT IMPLEMENTED
  assert(0);
}

