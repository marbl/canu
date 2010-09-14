#include "sim4polish.H"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "memory.h"

void
sim4polish::s4p_swapExons(u32bit a, u32bit b) {
  sim4polishExon  copyofa = _exons[a];

  _exons[a] = _exons[b];
  _exons[b] = copyofa;
}


//  Insert a single exon into the list at position a
void
sim4polish::s4p_insertExon(u32bit a, u32bit intronori, sim4polishExon *e) {
  sim4polish  p;

  p._numExons = 1;
  p._exons    = e;

  s4p_insertExons(a, intronori, &p);
}



//  Inserts all the exons in e into the list at position a.
void
sim4polish::s4p_insertExons(u32bit a, u32bit intronori, sim4polish *e) {
  sim4polishExon *ne = new sim4polishExon [_numExons + e->_numExons];

  //  Copy exons up to the insert point.

  for (u32bit i=0; i<a; i++) {
    ne[i] = _exons[i];
    _exons[i].s4p_clearExon();
  }

  //  Insert the new ones.  We don't own them, so can't assume anything about the alignment strings.

  for (u32bit i=0; i<e->_numExons; i++)
    ne[a+i].s4p_copyExon(e->_exons+i);

  //  Copy the rest.

  for (u32bit i=a; i<_numExons; i++) {
    ne[i+e->_numExons] = _exons[i];
    _exons[i].s4p_clearExon();
  }

  //  All done with the copy, get rid of the old stuff.  s4p_clearExon() above is critical here;
  //  without it we would delete the alignment strings.

  delete [] _exons;
  _exons = ne;

  _numExons += e->_numExons;

  //  We trust that the user has set the intron orientation in the new exon, and that 'intronori' is
  //  the correct orientation for the previous intron.
  //
  if (a > 0)
    _exons[a-1]._intronOrientation = intronori;
}
