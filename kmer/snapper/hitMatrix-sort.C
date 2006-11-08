#include "posix.H"
#include "hitMatrix.H"

//  Sort by dsPos

inline
void
adjustHeap_dsPos(diagonalLine *L, u32bit p, u32bit n) {
  u32bit  q = L[p]._qsPos;
  u32bit  d = L[p]._dsPos;
  u32bit  c = (p << 1) + 1;  //  let c be the left child of p

  while (c < n) {

    //  Find the larger of the two children
    //
    if ((c+1 < n) && (L[c]._dsPos < L[c+1]._dsPos))
      c++;

    //  Does the node in question fit here?
    //
    if (d >= L[c]._dsPos)
      break;

    //  Else, swap the parent and the child
    //
    L[p]._qsPos      = L[c]._qsPos;
    L[p]._dsPos      = L[c]._dsPos;

    //  Move down the tree
    //
    p = c;
    c = (p << 1) + 1;
  }

  L[p]._qsPos      = q;
  L[p]._dsPos      = d;
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
      u32bit  q  = _hits[i]._qsPos;
      u32bit  d  = _hits[i]._dsPos;

      _hits[i]._qsPos      = _hits[0]._qsPos;
      _hits[i]._dsPos      = _hits[0]._dsPos;

      _hits[0]._qsPos      = q;
      _hits[0]._dsPos      = d;
      
      adjustHeap_dsPos(_hits, 0, i);
    }
  }
}
