
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
 *  This file is derived from:
 *
 *    src/utility/AS_UTL_alloc.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-AUG-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef ARRAYS_H
#define ARRAYS_H

#include <algorithm>

using namespace std;


const uint32  resizeArray_doNothing = 0x00;
const uint32  resizeArray_copyData  = 0x01;
const uint32  resizeArray_clearNew  = 0x02;



template<typename TT, typename LL>
void
allocateArray(TT*& array, LL arrayMax, uint32 op=resizeArray_clearNew) {

  if (array != NULL)
    delete [] array;

  array = new TT [arrayMax];

  if (op == resizeArray_clearNew)
    memset(array, 0, sizeof(TT) * arrayMax);
}



template<typename TT>
TT *
duplicateString(TT const *fr) {

  if (fr == NULL)
    return(NULL);

  uint32  ln = strlen(fr);
  TT     *to = new TT [ln+1];

  memcpy(to, fr, sizeof(TT) * (ln+1));

  return(to);
}



template<typename TT, typename LL>
void
duplicateArray(TT*& to, LL &toLen, LL &toMax, TT *fr, LL frLen, LL frMax=0, bool forceAlloc=false) {

  if (fr == NULL)
    assert(frLen == 0);

  if ((toMax < frLen) || (forceAlloc)) {
    delete [] to;

    toMax = frLen;
    to    = new TT [toMax];
  }

  toLen = frLen;

  if (frLen > 0)
    memcpy(to, fr, sizeof(TT) * frLen);
}


//  Set the array size to 'newMax'.  No guards, the array will ALWAYS be reallocated.

template<typename TT, typename LL>
void
setArraySize(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, uint32 op=resizeArray_copyData) {

  arrayMax =     newMax;
  arrayLen = min(newMax, arrayLen);

  TT *copy = new TT [arrayMax];

  if ((op & resizeArray_copyData) && (array != NULL) && (arrayLen > 0))
    memcpy(copy, array, sizeof(TT) * arrayLen);

  delete [] array;
  array = copy;

  if ((op & resizeArray_clearNew) && (arrayMax > arrayLen))
    memset(array + arrayLen, 0, sizeof(TT) * (arrayMax - arrayLen));
}




//  Ensure that there is enough space to hold one more element in the array.
//  Increase the array by 'moreSpace' if needed.

template<typename TT, typename LL>
void
increaseArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 moreSpace) {
  uint64  newMax = arrayMax + moreSpace;

  if (arrayLen < arrayMax)
    return;

  setArraySize(array, arrayLen, arrayMax, newMax, resizeArray_copyData);
}


template<typename T1, typename T2, typename LL>
void
increaseArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, uint64 moreSpace) {
  uint64  newMax = arrayMax + moreSpace;

  if (arrayLen < arrayMax)
    return;

  setArraySize(array1, arrayLen, arrayMax, newMax, resizeArray_copyData);
  setArraySize(array2, arrayLen, arrayMax, newMax, resizeArray_copyData);
}


//  Resize the array so that it is at least as big as new max.  Do nothing
//  if the array is big enough already.

template<typename TT, typename LL>
void
resizeArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, uint32 op=resizeArray_copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array, arrayLen, arrayMax, newMax, op);
}


template<typename T1, typename T2, typename LL>
void
resizeArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, LL newMax, uint32 op=resizeArray_copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array1, arrayLen, arrayMax, newMax, op);
  setArraySize(array2, arrayLen, arrayMax, newMax, op);
}


#endif  //  ARRAYS_H
