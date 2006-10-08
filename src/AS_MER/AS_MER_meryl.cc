
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "AS_MER_stream.h"

static cds_int64 mersToCount = 999999999999LL;

////////////////////////////////////////////////////////////////////////////////
//
//  These are stolen from libbri/alphabet.C
//
unsigned char   compressSymbol[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   validSymbol[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   decompressSymbol[256] = { 65, 67, 71, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   complementSymbol[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84, 86, 71, 72, 0, 0, 67, 68, 0, 0, 77, 0, 75, 78, 0, 0, 0, 89, 87, 65, 65, 66, 83, 0, 82, 0, 0, 0, 0, 0, 0, 0, 116, 118, 103, 104, 0, 0, 99, 100, 0, 0, 109, 0, 107, 110, 0, 0, 0, 121, 119, 97, 97, 98, 115, 0, 114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   validCompressedSymbol[256] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };



////////////////////////////////////////////////////////////////////////////////
//
//  These are stolen from libbri/bit-packing.H
//    setDecodedValue
//    getDecodedValue
//    preDecrementDecodedValue
//
inline
cds_uint64
getDecodedValue(cds_uint64 *ptr,
                cds_uint64  pos,
                cds_uint64  siz) {
  cds_uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  cds_uint64 bit = (pos     ) & 0x000000000000003fllu;
  cds_uint64 b1  = 64 - bit;
  cds_uint64 b2  = siz - b1;  //  Only used if siz > b1
  cds_uint64 ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);
  } else {
    ret  = (ptr[wrd] & CDS_UINT64_MASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & CDS_UINT64_MASK(b2);
  }

  ret &= CDS_UINT64_MASK(siz);

  return(ret);
}

inline
void
setDecodedValue(cds_uint64 *ptr,
                cds_uint64  pos,
                cds_uint64  siz,
                cds_uint64  val) {
  cds_uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  cds_uint64 bit = (pos     ) & 0x000000000000003fllu;
  cds_uint64 b1  = 64 - bit;
  cds_uint64 b2  = siz - b1;  //  Only used if siz > b1

  val &= CDS_UINT64_MASK(siz);

  if (b1 >= siz) {
    ptr[wrd] &= ~( CDS_UINT64_MASK(siz) << (b1-siz) );
    ptr[wrd] |= val << (b1-siz);
  } else {
    ptr[wrd] &= ~CDS_UINT64_MASK(b1);
    ptr[wrd] |= (val & (CDS_UINT64_MASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(CDS_UINT64_MASK(b2) << (64-b2));
    ptr[wrd+1] |= (val & (CDS_UINT64_MASK(b2))) << (64-b2);
  }
}

inline
cds_uint64
preDecrementDecodedValue(cds_uint64 *ptr,
                         cds_uint64  pos,
                         cds_uint64  siz) {
  cds_uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  cds_uint64 bit = (pos     ) & 0x000000000000003fllu;
  cds_uint64 b1  = 64 - bit;
  cds_uint64 b2  = siz - b1;  //  Only used if siz > b1
  cds_uint64 ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= CDS_UINT64_MASK(siz);

    ptr[wrd] &= ~( CDS_UINT64_MASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & CDS_UINT64_MASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & CDS_UINT64_MASK(b2);

    ret--;
    ret &= CDS_UINT64_MASK(siz);

    ptr[wrd] &= ~CDS_UINT64_MASK(b1);
    ptr[wrd] |= (ret & (CDS_UINT64_MASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(CDS_UINT64_MASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (CDS_UINT64_MASK(b2))) << (64-b2);
  }

  return(ret);
}



//  Parameters for the mer counter, and a nice utility function.
//
class mcDescription {
public:
  cds_uint32      _merSizeInBases;
  cds_uint32      _merSizeInBits;

  cds_uint32      _tableSizeInBits;
  cds_uint64      _tableSizeInEntries;

  cds_uint32      _chckBits;
  cds_uint64      _chckMask;

  cds_uint32      _hashWidth;
  cds_uint64      _hashMask;

  cds_uint64      _actualNumberOfMers;

  cds_uint64 HASH(cds_uint64 a) {
    return((a >> _chckBits) & _hashMask);
  }
};

static mcDescription   mcd;

cds_uint64     *_chck;
cds_uint64     *_hash;

//  Stuff for sorting the list of mers.
//
typedef cds_uint64 heapbit;




const char *usage =
"usage: %s [options]\n"
"\n"
"Given a sequence file (-s) and lots of parameters, compute\n"
"the mer-count tables.  By default, both strands are processed.\n"
"\n"
"        -m #                    (size of a mer)\n"
"        -s /path/to/fragstore   (a fragstore)\n"
"        -n #                    (threshold to output)\n"
"        -N #                    (number of mers to examine)\n"
"        -K #                    (read only every Kth ffragment)\n"
"        -o tblprefix            (output table prefix)\n";





#if 0
cds_uint32
estimateTableSize(cds_uint64 numMers, cds_uint32 merSize) {

  //  How many bits do we need in the hash table to store all the mers?
  //
  cds_uint32 hBits = 1;
  while ((numMers+1) > (CDS_UINT64_ONE << hBits))
    hBits++;


  //  Compute the table size that results in the smalest memory footprint.
  //  Footprint is measured in megabytes.
  //
  cds_uint64 one  = CDS_UINT64_ONE;
  cds_uint64 cmin = 1048576;    //  Assume that using a terabyte is big.
  cds_uint64 t    = 0;

  for (cds_uint64 h=16; h<=32 && h<2*merSize; h++) {
    cds_uint64  c     = 2 * merSize - h;
    cds_uint64  hSize = ((one << h) * hBits) >> 23;
    cds_uint64  cSize = (numMers * c) >> 23;

    if (cmin > hSize+cSize) {
      cmin = hSize+cSize;
      t    = h;
    }
  }

  fprintf(stderr, "Estimated table size to be %d mers -> t = %d\n", numMers, t);

  return(t);
}
#endif










////////////////////////////////////////////////////////////////////////////////
//
//  BUILD
//
//

void
createHashTable(char *inputFile, cds_uint32 skipNum) {
  cds_uint64  mer;

  cds_uint32 *_ctbl = new cds_uint32 [ mcd._tableSizeInEntries ];
  for (cds_uint32 i=mcd._tableSizeInEntries; i--; )
    _ctbl[i] = 0;

  //
  //  XXX:  Huge stalls accessing _ctbl.
  //
  //  I think we can get by with 24 bits of count here; the fragstore
  //  seems to have only 9 million mers in the first bucket.  Should
  //  probably check for overflow, though.
  //  
  merStream          M(mcd._merSizeInBases, inputFile, skipNum);
  cds_int64 nummers=0;
  while (M.nextMer()&&(nummers++)<mersToCount) {
    mer = M.theFMer();
    if (mer > M.theRMer())
      mer = M.theRMer();

    _ctbl[ mcd.HASH(mer) ]++;



#if 0
    char theMerString[33];
    for (cds_uint32 i=0; i<mcd._merSizeInBases; i++)
      theMerString[mcd._merSizeInBases-i-1] = decompressSymbol[(M.theFMer() >> (2*i)) & 0x03];
    theMerString[mcd._merSizeInBases] = 0;
    fprintf(stdout, "%s\n", theMerString);
#endif

    mcd._actualNumberOfMers++;
  }


  //
  //  Allocate a PACKED array for the hash table.  This needs to be
  //  packed only for mcd._actualNumberOfMers > 4 billion, really.
  //

  //  Determine how many bits we need to hold the value
  //  mcd._actualNumberOfMers.....then....
  //
  //  This is mcd._actualNumberOfMers+1 because we need to store the
  //  first position after the last mer.  That is, if there are two
  //  mers, we will store that the first mer is at position 0, the
  //  second mer is at position 1, and the end of the second mer is at
  //  position 2.
  //
  mcd._hashWidth  = 1;
  while ((mcd._actualNumberOfMers+1) > (CDS_UINT64_ONE << mcd._hashWidth))
    mcd._hashWidth++;

  //  ....allocate a hash table that is that many bits wide.
  //
  _hash = new cds_uint64 [(mcd._tableSizeInEntries+1) * mcd._hashWidth / 64 + 2];

  //
  //  Create the hash index using the counts.  The hash points
  //  to the end of the bucket; when we add a word, we move the
  //  hash bucket pointer down one.
  //
  //  When done, we can deallocate the counting table.
  //
  cds_uint64 i=0;
  cds_uint64 j=0;
  cds_uint64 c=0;

  while (i < mcd._tableSizeInEntries) {
    c += _ctbl[i++];
    setDecodedValue(_hash, j, mcd._hashWidth, c);
    j += mcd._hashWidth;
  }

  //  Add the location of the end of the table.  This is not
  //  modified when adding words, but is used to determine
  //  the size of the last bucket.
  //
  setDecodedValue(_hash, j, mcd._hashWidth, c);

  delete [] _ctbl;
}


void
verifyHashTable(void) {
  cds_uint64 i=0, j=0, c=0, d=0;

  fprintf(stderr, "    Verifying hash table.\n");

  while (i <= mcd._tableSizeInEntries) {
    c = getDecodedValue(_hash, j, mcd._hashWidth);

    if (c < d)
      fprintf(stderr, "ERROR:  Table[%lu] out of order.\n", i);

    d = c;

    j += mcd._hashWidth;
    i++;
  }

  fprintf(stderr, "    Verify finished.\n");
}


//  0.52 Mb per second on viking5
void
fillCheckTable(char *inputFile, cds_uint32 skipNum) {
  cds_uint64  mer, b, c;

  //  Allocate space for mcd._actualNumberOfMers mers in the _chck array.
  //  This doesn't need to be cleared.
  //
  _chck = new cds_uint64 [mcd._actualNumberOfMers * mcd._chckBits / 64 + 1];

  merStream   M(mcd._merSizeInBases, inputFile,skipNum);
  cds_int64 nummers=0;

  while (M.nextMer()&&(nummers++)<mersToCount) {
    mer = M.theFMer();
    if (mer > M.theRMer())
      mer = M.theRMer();

    b = mcd.HASH(mer) * mcd._hashWidth;
    c = preDecrementDecodedValue(_hash, b, mcd._hashWidth) * mcd._chckBits;

    setDecodedValue(_chck, c, mcd._chckBits, mer & mcd._chckMask);
  }
}






////////////////////////////////////////////////////////////////////////////////
//
//  OUTPUT
//
void
adjustHeap(heapbit *M, cds_int64 i, cds_int64 n) {
  heapbit   m = M[i];
  cds_int64    j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && M[j] < M[j+1])
      j++;                   //  j is the larger child

    if (m >= M[j])           //  a position for M[i] has been found
      break;

    M[(j-1)/2] = M[j];       //  Move larger child up a level

    j = (j << 1) + 1;
  }

  M[(j-1)/2] = m;
}



void
sortAndOutput(char   *outfilename,
              cds_uint32  targetCount) {
  cds_uint64 m     = CDS_UINT64_ONE << mcd._tableSizeInBits;
  cds_uint32 count = 0;
  cds_uint32 items = 0;

  heapbit *_sortedList    = 0L;
  cds_uint32   _sortedListMax = 0;
  cds_uint32   _sortedListLen = 0;

  unsigned char  theMerString[33];
  cds_uint64         mer;

  //  Open the output file
  //
  FILE *outFile = fopen(outfilename, "w");


  //  For each bucket, sort it.  The output is done
  //  in the sort.
  //
  for (cds_uint64 B=0, b=0; b<m; b++) {
    cds_uint64 st = getDecodedValue(_hash, B, mcd._hashWidth);
    B        += mcd._hashWidth;
    cds_uint64 ed = getDecodedValue(_hash, B, mcd._hashWidth);

    if (ed < st) {
      fprintf(stderr, "ERROR: Bucket %10lu starts at %10lu ends at %10lu\n", b, st, ed);
      fflush(stderr);
    }

    _sortedListLen = ed - st;

    count = 0;
    items = 0;

    if (_sortedListLen > 0) {

      //  Allocate more space, if we need to.
      //
      if (_sortedListLen > _sortedListMax) {
        delete [] _sortedList;
        _sortedList    = new heapbit [_sortedListLen + 1];
        _sortedListMax = _sortedListLen;
      }

      //  Unpack the check values
      //
      for (cds_uint64 i=st, J=st*mcd._chckBits; i<ed; i++, J += mcd._chckBits)
        _sortedList[i-st] = getDecodedValue(_chck, J, mcd._chckBits);

      //  Sort if there is more than one item
      //
      if (_sortedListLen > 1) {

        //  Create the heap of lines.
        //
        for (cds_int64 t=(_sortedListLen-2)/2; t>=0; t--)
          adjustHeap(_sortedList, t, _sortedListLen);

        //  Interchange the new maximum with the element at the end of the tree
        //
        for (cds_int64 t=_sortedListLen-1; t>0; t--) {
          heapbit          tv = _sortedList[t];
          _sortedList[t]      = _sortedList[0];
          _sortedList[0]      = tv;

          adjustHeap(_sortedList, 0, t);
        }
      }


      //  Scan the list of sorted mers, counting them.  Whenever we 
      //  know the count, output it.
      //
      count = 1;
      if (_sortedListLen > 0) {
        cds_uint32 t;
        for (t=1; t<_sortedListLen; t++) {
          if (_sortedList[t] != _sortedList[t-1]) {
            if (targetCount <= count) {
              mer = (b << mcd._chckBits) | _sortedList[t-1];
              for (cds_uint32 i=0; i<mcd._merSizeInBases; i++)
                theMerString[mcd._merSizeInBases-i-1] = decompressSymbol[(mer >> (2*i)) & 0x03];
              theMerString[mcd._merSizeInBases] = 0;
              fprintf(outFile, ">%d\n%s\n", count, theMerString);
            }
            count = 0;
          }

          count++;
        }

        //  Dump the last mer
        //
        if (targetCount <= count) {
          mer = (b << mcd._chckBits) | _sortedList[t-1];
          for (cds_uint32 i=0; i<mcd._merSizeInBases; i++)
            theMerString[mcd._merSizeInBases-i-1] = decompressSymbol[(mer >> (2*i)) & 0x03];
          theMerString[mcd._merSizeInBases] = 0;
          fprintf(outFile, ">%d\n%s\n", count, theMerString);
        }
      }
    }
  }

  fclose(outFile);
}






void
build(char   *fragStore,
      char   *outputFile,
      cds_uint32  merSize,
      cds_uint32  targetCount,
      cds_uint32 skipNum) {

  mcd._merSizeInBases      = merSize;
  mcd._merSizeInBits       = mcd._merSizeInBases << 1;

  //  XXX:  Hardcoded; need to estimate based on the fragstore
  //  estimateTableSize(fragStore, merSize),
  //
  mcd._tableSizeInBits     = 26;
  mcd._tableSizeInEntries  = 1 << mcd._tableSizeInBits;

  mcd._chckBits            = mcd._merSizeInBits - mcd._tableSizeInBits;
  _chck                    = 0L;
  mcd._chckMask            = CDS_UINT64_MASK(mcd._chckBits);

  mcd._hashWidth           = 0;
  _hash                    = 0L;
  mcd._hashMask            = CDS_UINT64_MASK(mcd._tableSizeInBits);  //  unused?

  mcd._actualNumberOfMers  = 0;

  createHashTable(fragStore,skipNum);
  verifyHashTable();
  fillCheckTable(fragStore,skipNum);
  sortAndOutput(outputFile, targetCount);

  delete [] _chck;
  delete [] _hash;
}












int
main(int argc, char **argv) {
  cds_uint32            merSize          = 20;
  char             *fragStore        = 0L;
  char             *outputFile       = 0L;
  cds_uint64            minimumCount     = 0;
  cds_uint32            skipNum = 1;

  if (argc == 1) {
    fprintf(stderr, usage, argv[0], argv[0]);
    exit(1);
  }

  for (int arg=1; arg < argc; arg++) {
    if (argv[arg][0] != '-') {
      fprintf(stderr, "Not an option: '%s'.\n", argv[arg]);
      exit(1);
    } else {
      switch (argv[arg][1]) {
        case 'm':
          arg++;
          merSize = atoi(argv[arg]);
          break;
        case 's':
          arg++;
          fragStore = argv[arg];
          break;
        case 'n':
          arg++;
          minimumCount = STR_TO_UINT64(argv[arg], NULL, 10);
          break;
        case 'K':
	  arg++;
	  skipNum = atoi(argv[arg]);
	  break;
        case 'N':
	  arg++;
          mersToCount = STR_TO_UINT64(argv[arg], NULL, 10);
 	  if(mersToCount<=0){
 	    fprintf(stderr,"Trouble getting number of mers to count from %s (option -N)\n",
 		    argv[arg]);
 	    exit(1);
 	  }
	  break;
        case 'o':
          arg++;
          outputFile = argv[arg];
          break;
        default:
          fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
          break;
      }
    }
  }

  if (fragStore == 0L) {
    fprintf(stderr, "ERROR - no fragstore loaction specified.\n");
    exit(1);
  }

  if (outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }

  build(fragStore,
        outputFile,
        merSize,
        minimumCount,
	skipNum);
}
