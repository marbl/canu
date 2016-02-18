
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
 *    Brian P. Walenz from 2005-AUG-24 to 2013-AUG-01
 *      are Copyright 2005-2008,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-11 to 2015-JAN-31
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-08
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_ALLOC_H
#define AS_UTL_ALLOC_H


template<typename TT, typename LL>
void
duplicateArray(TT*& to, LL &toLen, LL &toMax, TT *fr, LL frLen, LL UNUSED(frMax), bool forceAlloc=false) {

  if ((toMax < frLen) || (forceAlloc)) {
    delete [] to;

    toLen = frLen;
    toMax = toLen;

    to = new TT [toMax];
  }

  memcpy(to, fr, sizeof(TT) * toLen);
}

const uint32  resizeArray_doNothing = 0x00;
const uint32  resizeArray_copyData  = 0x01;
const uint32  resizeArray_clearNew  = 0x02;


template<typename TT, typename LL>
void
resizeArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, uint32 op=resizeArray_copyData) {

  if (newMax <= arrayMax)
    return;

  arrayMax = newMax;

  TT *copy = new TT [arrayMax];

  if (op & resizeArray_copyData)
    memcpy(copy, array, sizeof(TT) * arrayLen);

  delete [] array;
  array = copy;

  if (op & resizeArray_clearNew)
    memset(array + arrayLen, 0, sizeof(TT) * (newMax - arrayLen));
}


template<typename T1, typename T2, typename LL>
void
resizeArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, LL newMax, uint32 op=resizeArray_copyData) {

  if (newMax <= arrayMax)
    return;

  arrayMax = newMax;

  T1 *copy1 = new T1 [arrayMax];
  T2 *copy2 = new T2 [arrayMax];

  if (op & resizeArray_copyData) {
    memcpy(copy1, array1, sizeof(T1) * arrayLen);
    memcpy(copy2, array2, sizeof(T2) * arrayLen);
  }

  delete [] array1;
  delete [] array2;

  array1 = copy1;
  array2 = copy2;

  if (op & resizeArray_clearNew) {
    memset(array1 + arrayLen, 0, sizeof(T1) * (newMax - arrayLen));
    memset(array2 + arrayLen, 0, sizeof(T2) * (newMax - arrayLen));
  }
}


template<typename TT, typename LL>
void
increaseArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 increment) {

  if (arrayLen + increment <= arrayMax)
    return;

  LL newMax = arrayMax;

  while (newMax < arrayLen + increment)
    newMax *= 2;

  resizeArray(array, arrayLen, arrayMax, newMax, resizeArray_copyData);
}


template<typename T1, typename T2, typename LL>
void
increaseArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, uint64 increment) {

  if (arrayLen + increment <= arrayMax)
    return;

  LL newMax = arrayMax;

  while (newMax < arrayLen + increment)
    newMax *= 2;

  resizeArrayPair(array1, array2, arrayLen, arrayMax, newMax, resizeArray_copyData);
}


#endif // AS_UTL_ALLOC_H
