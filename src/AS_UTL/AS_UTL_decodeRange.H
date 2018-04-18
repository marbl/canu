
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
 *    Brian P. Walenz beginning on 2018-APR-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
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
char *
AS_UTL_decodeRange(char *range, T &lo, T &hi) {
  char    *ap = range;

  strtonumber(lo, ap, &ap);     //  Grab the first number.

  hi = lo;                      //  Set the second to that.

  if (*ap == '-') {             //  If this is a range,
    ap++;                       //  grab the second number
    strtonumber(hi, ap, &ap);
  }

  if (*ap == ',')               //  If the next letter continues
    return(ap + 1);             //  move past that and return.

  if (*ap == 0)                 //  If the next letter is the end
    return(NULL);               //  of the string, return NULL.

  //  Otherwise, we can't decode this range.

  fprintf(stderr, "ERROR: invalid range '%s'\n", range);
  exit(1);

  return(NULL);
}



template<typename T>
void
AS_UTL_decodeRange(char *range, vector<T> &bgn, vector<T> &end) {
  char    *ap = range;
  T        av = 0;
  T        bv = 0;

  while ((ap != NULL) && (*ap != 0)) {
    range = AS_UTL_decodeRange(range, av, bv);

    bgn.push_back(av);
    end.push_back(bv);
  }
}



template<typename T>
void
AS_UTL_decodeRange(char *range, set<T> &ranges) {
  char    *ap = range;
  T        av = 0;
  T        bv = 0;

  while ((ap != NULL) && (*ap != 0)) {
    range = AS_UTL_decodeRange(range, av, bv);

    for (T xx=av; xx<=bv; xx++)
      ranges.insert(xx);
  }
}


#endif  //  AS_UTL_DECODERANGE_H
