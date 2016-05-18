
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
 *    src/bogart/AS_BAT_Datatypes.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-APR-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_FRAGMENT_INFO
#define INCLUDE_AS_BAT_FRAGMENT_INFO

#include "AS_global.H"
#include "ovStore.H"
#include "gkStore.H"

#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>



class FragmentInfo {
public:
  FragmentInfo(gkStore *gkp, const char *prefix, uint32 minReadLen);
  ~FragmentInfo();

  uint64  memoryUsage(void) {
    return((3 * sizeof(uint32) * _numFragments) +
           (2 * sizeof(double) * _numLibraries) +
           (2 * sizeof(uint32) * _numLibraries));
  };

  uint32  numFragments(void) { return(_numFragments); };
  uint32  numLibraries(void) { return(_numLibraries); };

  uint32  fragmentLength(uint32 iid) { return(_fragLength[iid]); };
  uint32  libraryIID(uint32 iid)     { return(_libIID[iid]);  };

  uint32  overlapLength(uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang) {
    int32  alen = fragmentLength(a_iid);
    int32  blen = fragmentLength(b_iid);
    int32  aovl = 0;
    int32  bovl = 0;

    assert(alen > 0);
    assert(blen > 0);

    if (a_hang < 0) {
      //  b_hang < 0      ?     ----------  :     ----
      //                  ?  ----------     :  ----------
      //
      aovl = (b_hang < 0) ? (alen + b_hang) : (alen);
      bovl = (b_hang < 0) ? (blen + a_hang) : (blen + a_hang - b_hang);
    } else {
      //  b_hang < 0      ?  ----------              :  ----------
      //                  ?     ----                 :     ----------
      //
      aovl = (b_hang < 0) ? (alen - a_hang + b_hang) : (alen - a_hang);
      bovl = (b_hang < 0) ? (blen)                   : (blen - b_hang);
    }

    if ((aovl <= 0) || (bovl <= 0) || (aovl > alen) || (bovl > blen)) {
      fprintf(stderr, "WARNING: bogus overlap found for A="F_U32" B="F_U32"\n", a_iid, b_iid);
      fprintf(stderr, "WARNING:                     A len="F_S32" hang="F_S32" ovl="F_S32"\n", alen, a_hang, aovl);
      fprintf(stderr, "WARNING:                     B len="F_S32" hang="F_S32" ovl="F_S32"\n", blen, b_hang, bovl);
    }

    if (aovl < 0)     aovl = 0;
    if (bovl < 0)     bovl = 0;

    if (aovl > alen)  aovl = alen;
    if (bovl > blen)  bovl = blen;

    assert(aovl > 0);
    assert(bovl > 0);
    assert(aovl <= alen);
    assert(bovl <= blen);

    //  AVE does not work.      return((uint32)((aovl, bovl)/2));
    //  MAX does not work.      return((uint32)MAX(aovl, bovl));

    return(aovl);
  };

private:
  void      save(const char *prefix);
  bool      load(const char *prefix);

  uint32   _numFragments;
  uint32   _numLibraries;

  uint32  *_fragLength;
  uint32  *_libIID;
};



extern FragmentInfo     *FI;

#endif  //  INCLUDE_AS_BAT_FRAGMENT_INFO
