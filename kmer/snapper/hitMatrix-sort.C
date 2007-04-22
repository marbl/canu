#include "snapper2.H"

//  Sort by dsPos

inline
void
adjustHeap_dsPos(diagonalLine *L, u32bit p, u32bit n) {
  u64bit  v = L[p].all;
  u64bit  d = L[p].val.dPos;
  u32bit  c = (p << 1) + 1;  //  let c be the left child of p

  while (c < n) {

    //  Find the larger of the two children
    //
    if ((c+1 < n) && (L[c].val.dPos < L[c+1].val.dPos))
      c++;

    //  Does the node in question fit here?
    //
    if (d >= L[c].val.dPos)
      break;

    //  Else, swap the parent and the child
    //
    L[p].all = L[c].all;

    //  Move down the tree
    //
    p = c;
    c = (p << 1) + 1;
  }

  L[p].all = v;
}

void
hitMatrix::sort_dsPos(void) {

  if (_hitsLen > 1) {

    //  Create the heap of lines.
    //
    for (u32bit i=_hitsLen/2; i--; )
      adjustHeap_dsPos(_hits, i, _hitsLen);

    //  Interchange the new maximum with the element at the end of the tree
    //
    for (u32bit i=_hitsLen-1; i>0; i--) {
      u64bit v       = _hits[i].all;
      _hits[i].all   = _hits[0].all;
      _hits[0].all   = v;
      
      adjustHeap_dsPos(_hits, 0, i);
    }
  }
}
