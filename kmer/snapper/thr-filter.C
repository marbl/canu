#include "thr.H"



u32bit
configureFilter(double L,
                double H,
                double V,
                aHit  *theHits,
                u32bit theHitsLen) {

  //  Find the highest and lowest quality hit
  //
  u32bit  hiQ = theHits[0]._covered;
  u32bit  loQ = hiQ;

  for (u32bit i=0; i<theHitsLen; i++) {
    if (hiQ < theHits[i]._covered)
      hiQ = theHits[i]._covered;
    if (loQ > theHits[i]._covered)
      loQ = theHits[i]._covered;
  }

  //  _numMers is not the same as the number covered, so we should
  //  ensure that h is in range.
  //
  //  Note: _numMers is constant for all hits, so we can use any of them
  //
  double h = (double)(hiQ - loQ) / (double)theHits[0]._numMers;
  if (h > 1.0)
    h = 1.0;

  double p = 0.0;
  if      (h <= L)  p = 1.0;
  else if (h >= H)  p = V;
  else              p = 1.0 - (1.0 - V) * (h - L) / (H - L);

  if (p > 1.0) {
    fprintf(stderr, "error in p; p=%f > 1.0!  h=%f (L=%f H=%f V=%f)\n", p, h, L, H, V);
    p = 1.0;
  }

  if (V - p > 1e-10) {
    fprintf(stderr, "error in p; p=%f < V!  h=%f (L=%f H=%f V=%f)\n", p, h, L, H, V);
    p = V;
  }

  //  Any thing at or above cutL is good, and we should polish it.
  //  Anything below is junk, and we should ignore it.
  //
  return((u32bit)floor(hiQ - p * h * theHits[0]._numMers));
}



u32bit
doFilter(searcherState       *state,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         logMsg              *theLog) {

  if (theHitsLen == 0)
    return(0);

  u32bit numF = 0;
  u32bit cutL = configureFilter(config._Lo,
                                config._Hi,
                                config._Va, theHits, theHitsLen);

  for (u32bit i=0; i < theHitsLen; i++) {

    //  If the coverage of the hit is more than the minimum, mark the
    //  hit as polishable.  Unless the hit was discarded.
    //
    if (!(theHits[i]._status & AHIT_DISCARDED) &&
        (theHits[i]._covered >= cutL)) {
      numF++;
      theHits[i]._status |= AHIT_POLISHABLE;
    }
  }

  return(numF);
}
