
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
 *    kmer/libutil/u32bitList.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2006-JUN-13 to 2014-APR-11
 *      are Copyright 2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef UINT32LIST_H
#define UINT32LIST_H

#include <string.h>

//  A very simple integer list.  Hopefully lighter weight than a
//  vector.

//  It might be useful to extend this to have 'undef' values,
//  and to allow shift(), pop().

class uint32List {
public:
  uint32List(uint32 max=16) {
    _len = 0;
    _max = max;
    _lst = new uint32 [_max];
  };
  ~uint32List() {
    delete [] _lst;
  };

private:
  void     resize(uint32 idx) {
    if (idx >= _max) {
      _max *= 2;
      uint32 *L = new uint32 [_max];
      memcpy(L, _lst, sizeof(uint32) * _len);
      delete [] _lst;
      _lst = L;
    }
    if (idx >= _len)
      _len = idx + 1;
  }

public:
  uint32 &operator[](uint32 idx) {
    resize(idx);
    return(_lst[idx]);
  }

  void    push(uint32 val) {
    resize(_len);
    _lst[_len++] = val;
  }

  uint32  length(void) {
    return(_len);
  };

  void    clear(void) {
    _len = 0;
  }

private:
  uint32  _len;
  uint32  _max;
  uint32 *_lst;
};


#endif  //   UINT32LIST_H
