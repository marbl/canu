#ifndef MER_COVERING_H
#define MER_COVERING_H

//  This is an interval list, where the intervals are built using
//  fixed size pieces.
//
//  It's designed to accept pieces in roughly sorted order.
//
//  Intervals are stored c-style.
//

#include <stdio.h>
#include <stdlib.h>

class merCovering {
private:
  class interval {
  public:
    uint32      _lo;
    uint32      _hi;
    interval   *_next;

    interval(uint32 lo, uint32 hi, interval *n) {
      _lo   = lo;
      _hi   = hi;
      _next = n;
    }
  };

  interval    *_intervals;
  uint32       _width;
  uint32       _pieces;
#ifdef TEST_MERCOVERING
  uint32       _test[TEST_SIZE];
#endif

public:
  merCovering(uint32 w) {
    _intervals = 0L;
    _width     = w;
    _pieces    = 0;

#ifdef TEST_MERCOVERING
    for (uint32 i=0; i<TEST_SIZE; i++)
      _test[i] = 0;
#endif
  };

  ~merCovering() {
    clear();
  };

  void          clear(void) {
    interval *i = _intervals;

    while (i) {
      _intervals = i->_next;
      delete i;
      i = _intervals;
    }

    _intervals = 0L;
    _pieces    = 0;
  };

  uint32        sumOfLengths(void) {
    uint32 s=0;

    for (interval *i=_intervals; i; i = i->_next)
      s += i->_hi - i->_lo;

    return(s);
  };

  uint32        numberOfPieces(void) {
    return(_pieces);
  };

  void          addMer(uint32 lo) {
    _pieces++;

    uint32 hi = lo + _width;

    interval  *c;

#ifdef TEST_MERCOVERING
    for (uint32 i=lo; i<hi; i++)
      _test[i] = 1;
#endif

    //  Case:  No existing intervals, or the new interval extends
    //  the low range of the first interval
    //
    if ((_intervals == 0L) ||
        (hi < _intervals->_lo)) {
      _intervals = new interval(lo, hi, _intervals);
      return;
    }

    c = _intervals;

    while (c) {

      //  Case:  New interval is completely contained in the current interval.
      //
      if ((c->_lo <= lo) && (hi <= c->_hi))
        return;

      //  Case:  New interval overlaps the low end of the current interval,
      //  or is completely contained in an existing interval.
      //
      if ((lo <= c->_lo) && (hi <= c->_hi)) {
        c->_lo = lo;
        return;
      }

      if (c->_next) {

        //  Case: New interval overlaps the high end of the current interval...
        //
        if (lo <= c->_hi) {

          if (hi < c->_next->_lo) {
            //  but does not intersect the next interval.
            //
            c->_hi = hi;
            return;
          } else {
            //  and does intersect the next interval.
            //
            interval *p = c->_next;

            c->_hi   = c->_next->_hi;
            c->_next = c->_next->_next;

            delete p;
            return;
          }
        } else {
          //  Case: New interval is between two existing intervals
          //
          //  (lo > c->_hi) is given
          //
          if (hi < c->_next->_lo) {
            c->_next = new interval(lo, hi, c->_next);
            return;
          }
        }
      } else {
        //  Case:  New interval overlaps the high end of the current interval
        //
        if (lo <= c->_hi) {
            c->_hi = hi;
            return;
        } else {
          //  Otherwise, we just fell off the end of all intervals.
          //  Add one at the end.
          //
          c->_next = new interval(lo, hi,0L);
          return;
        }
      }

      c = c->_next;
    }

#ifdef TEST_MERCOVERING
    fprintf(stderr, "ERROR IN addInterval!\n");
#endif
  };

#ifdef TEST_MERCOVERING
  void         test(void) {
    for (uint32 i=0; i<TEST_SIZE; i++) {
      if (_test[i])
        _test[i] = 2;
    }
    for (interval *z=_intervals; z; z = z->_next) {
      for (uint32 i=z->_lo; i<z->_hi; i++) {
        if (_test[i] == 0) {
          fprintf(stderr, "INTERVAL CONTAINS SOMETHING NOT IN ARRAY! (%d)\n", i);
          exit(1);
        }
        if (_test[i] == 1) {
          fprintf(stderr, "INTERVAL HIT SOMETHING TWICE! (%d)\n", i);
          exit(1);
        }
        _test[i] = 1;
      }
    }
    for (uint32 i=0; i<TEST_SIZE; i++) {
      if (_test[i] == 2) {
        fprintf(stderr, "ARRAY CONTAINED SOMETHING NOT IN INTERVAL! (%d)\n", i);
        exit(1);
      }
    }
  };
#endif



  //  Incorporates the intervals in B into our list.
  //
  void         merge(merCovering *I = 0L) {
    interval *A, *B, *N, *L;

    if (I == 0L)
      return;

    A = _intervals;
    B = I->_intervals;
    N = 0L;
    L = 0L;

    while (A || B) {
      uint32 lo = 0;
      uint32 hi = 0;

      //  if either list is zero, we can just zip down the other list
      //  and add things.
      //
      if (!B) {
        while (A) {
          L->_next = new interval(A->_lo, A->_hi, 0L);
          L = L->_next;
          A = A->_next;
        }
      }

      if (!A) {
        while (B) {
          L->_next = new interval(B->_lo, B->_hi, 0L);
          L = L->_next;
          B = B->_next;
        }
      }

      if (A && B) {
        if (A->_lo == B->_lo) {
          //  A and B start at the same position
          //
          lo = A->_lo;
          hi = A->_hi;
          if (hi < B->_hi)
            hi = B->_hi;

          A = A->_next;
          B = B->_next;
        } else {
          //  A and B start at different positions.  Pick the first one.
          //
          if (A->_lo < B->_lo) {
            lo = A->_lo;
            hi = A->_hi;
            A = A->_next;
          } else {
            lo = B->_lo;
            hi = B->_hi;
            B = B->_next;
          }
        }

        //  We have an initial interval.  Add more stuff, while there
        //  are overlaps.

        bool modified = true;

        while ((A || B) && (modified)) {
          modified = false;

          if ((A) && (hi >= A->_lo)) {
            if (hi < A->_hi)
              hi = A->_hi;
            A = A->_next;
            modified = true;
          }

          if ((B) && (hi >= B->_lo)) {
            if (hi < B->_hi)
              hi = B->_hi;
            B = B->_next;
            modified = true;
          }
        }

        //  OK, got the new interval.  Save it.
        //
        if (N) {
          L->_next = new interval(lo, hi, 0L);
          L = L->_next;
        } else {
          N = L = new interval(lo, hi, 0L);
        }
      }
    }

    //  Save the number of mers in both intervals
    //
    uint32  p = _pieces + I->_pieces;

    clear();

    _intervals = N;
    _pieces    = p;
  }




#ifdef TEST_MERCOVERING
  void         dump(void) {
    for (interval *i=_intervals; i; i = i->_next)
      fprintf(stderr, "%5d-%5d ", i->_lo, i->_hi);
    fprintf(stderr, "\n");
  };

  void         compare(merCovering *B) {
    interval *i = _intervals;
    interval *j = B->_intervals;

    if (_pieces != B->_pieces) {
      fprintf(stderr, "Pieces differ (this=%d that=%d).\n", _pieces, B->_pieces);
      exit(1);
    }

    while (i && j) {
      if ((i->_lo != j->_lo) || (i->_hi != j->_hi)) {
        fprintf(stderr, "ERROR!\n");
        exit(1);
      }

      i = i->_next;
      j = j->_next;
    }

    if (i) {
      fprintf(stderr, "ERROR (i still exists)!\n");
      exit(1);
    }

    if (j) {
      fprintf(stderr, "ERROR (i still exists)!\n");
      exit(1);
    }
  };
#endif

};

#endif  //  MERCOVERING_H



