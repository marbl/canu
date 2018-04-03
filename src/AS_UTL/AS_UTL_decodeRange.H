
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-FEB-12 to 2013-OCT-11
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-MAY-28
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_DECODERANGE_H
#define AS_UTL_DECODERANGE_H

#include "AS_global.H"

#include <set>
#include <vector>

using namespace std;



template<typename T>
class AS_UTL_range {
public:
  AS_UTL_range() {
  };
  AS_UTL_range(T bgn, T end) {
  }

  bool  operator<(const AS_UTL_range<T> &that) const {
    if (_bgn < that._bgn)
      return(true);
    return(_end < that._end);
  };

private:
  T  _bgn;
  T  _end;
};



template<typename T>
void
AS_UTL_decodeRange(char *range, vector< AS_UTL_range<T> > &ranges) {
  char    *ap = range;
  T        av = 0;
  T        bv = 0;

  while (*ap != 0) {
    strtonumber(av, ap, &ap);

    if (*ap == ',') {
      ap++;
      ranges.insert(av);

    } else if (*ap == 0) {
      ranges.insert(av);

    } else if (*ap == '-') {
      ap++;
      strtonumber(bv, ap, &ap);

      for (T xx=av; xx<=bv; xx++)
        ranges.insert(xx);

      if (*ap == ',')
        ap++;

    } else if (*ap != 0) {
      fprintf(stderr, "ERROR: invalid range '%s'\n", range);
      exit(1);
    }
  }
}



template<typename T>
void
AS_UTL_decodeRange(char *range, set<T> &ranges) {
  char    *ap = range;
  T        av = 0;
  T        bv = 0;

  while (*ap != 0) {
    strtonumber(av, ap, &ap);

    if (*ap == ',') {
      ap++;
      ranges.insert(av);

    } else if (*ap == 0) {
      ranges.insert(av);

    } else if (*ap == '-') {
      ap++;
      strtonumber(bv, ap, &ap);

      for (T xx=av; xx<=bv; xx++)
        ranges.insert(xx);

      if (*ap == ',')
        ap++;

    } else if (*ap != 0) {
      fprintf(stderr, "ERROR: invalid range '%s'\n", range);
      exit(1);
    }
  }
}



template<typename T>
void
AS_UTL_decodeRange(char *range, T &lo, T &hi) {
  char    *ap = range;

  strtonumber(lo, ap, &ap);
  hi = lo;

  if (*ap == '-') {
    ap++;
    strtonumber(hi, ap, &ap);

  } else if (*ap != 0) {
    fprintf(stderr, "ERROR: invalid range '%s'\n", range);
    exit(1);
  }
}


#endif  //  AS_UTL_DECODERANGE_H
