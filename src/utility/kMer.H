
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
 *    kmer/libbio/kmer.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2005-MAY-19 to 2014-APR-11
 *      are Copyright 2005-2008,2010-2011,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received (LICENSE.txt) a copy of the GNU General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//  A 'simple' kMer datastructure.

#ifndef BIO_KMER_H
#define BIO_KMER_H

//  The maximum size of a mer.  You get 32 bases per word, so
//  KMER_WORDS=4 will get you up to a 128-mer.
//
#define KMER_WORDS  1

#include "AS_global.H"

#include "bitPackedFile.H"
#include "dnaAlphabets.H"

//  The kmer class is in two mutually exclusive implementations.  The first implementation
//  uses one word to store up to 32-mers.  The second implementation uses multiple words
//  to store arbitrarily large kmers.
//
//  Incomplete and out of date description of the methods:
//
//  //  Reverse all the words, reverse and complement the bases in
//  //  each word, then shift right to align the edge.
//  //
//  virtual kMerInterface  &reverseComplement(void) = 0;
//  virtual void            clear(void);
//
//  //  Construct a mer by shifting bases onto the end:
//  //    += shifts onto the right end
//  //    -= shifts onto the left end
//  //
//  virtual void   operator+=(uint64 x) = 0;
//  virtual void   operator-=(uint64 x) = 0;
//
//  //  used by merStream at least
//  //
//  virtual void     mask(bool) = 0;
//
//  //  Return the mer, as a 64-bit integer.  If the mer is more than
//  //  32-bases long, then the left-most (the earliest, the start, etc)
//  //  bases are used.
//  //
//  virtual operator uint64 () const = 0;
//
//  //  These are written/read in 5'endian, which isn't the most natural
//  //  implementation.  It's done this way to keep the sequence in
//  //  order (e.g., the merStreamFile).  Don't change the order.
//  //
//  //  On the otherhand, the implementation (of write anyway) is
//  //  basically the same as merToString().
//  //
//  //  Takes an optional number of BITS to write, pulled from the
//  //  END of the mer.
//  //
//  virtual void   writeToBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) const = 0;
//  virtual void   readFromBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) = 0;
//
//  //  Returns a sub-mer from either the start (left end) or the end
//  //  (right end) of the mer.  The sub-mer must be at most 64 bits
//  //  long.  Yes, BITS.
//  //
//  //  The start is difficult, because it can span multiple words.  The
//  //  end is always in the first word.
//  //
//  virtual uint64   startOfMer(uint32 bits) const = 0;
//  virtual uint64   endOfMer(uint32 bits) const = 0;
//
//  //  Set 'numbits' bits from (the end of) 'val' at bit position 'pos'
//  //  in the mer.  This is wildly low-level, but merylStreamReader
//  //  needs it.
//  //
//  //  The position is measured from the right end.
//  //    (0, 8, X) would copy the bits 7 to 0 of X to bits 7 to 0 of the mer.
//  //
//  //  Argh!  Can't use set/getDecodedValue because that is doing things in the wrong order.
//  //
//  //  Meryl
//  //
//  virtual uint64   getWord(uint32 wrd) const        = 0;  // { return(MERWORD(wrd)); };
//  virtual void     setWord(uint32 wrd, uint64 val)  = 0;  // { MERWORD(wrd) = val;   };
//
//  //  Show the mer as ascii
//  //
//  //  Doesn't print the last full word, if it's on the word boundary
//  //
//  //  We build the string right to left, print any partial word first,
//  //  then print whole words until we run out of words to print.
//  //
//  virtual char    *merToString(char *instr) const = 0;
//

#if KMER_WORDS == 1

#include "kMerTiny.H"
typedef kMerTiny kMer;

#else

#include "kMerHuge.H"
typedef kMerHuge kMer;

#endif


#undef DEBUGADDBASE
#undef DEBUGCOMP
#undef DEBUGSPACE


class kMerBuilder {
public:
  kMerBuilder(uint32 ms=0, uint32 cm=0, char *tm=0L);
  ~kMerBuilder();

  //  Clear all mer data, reset state to as just after construction.
  void    clear(bool clearMer=true);

  //  Returns true if we need another base to finish the mer.  This
  //  only occurs for compressed mers, if we are in a homopolymer run.
  //
private:
  bool addBaseContiguous(uint64 cf, uint64 cr);
  bool addBaseCompressed(uint64 cf, uint64 cr);
  bool addBaseSpaced(uint64 cf, uint64 cr);
  bool addBaseCompressedSpaced(uint64 cf, uint64 cr);

public:
  bool    addBase(char ch) {
    uint64  cf = alphabet.letterToBits(ch);
    uint64  cr = alphabet.letterToBits(alphabet.complementSymbol(ch));

#ifdef DEBUGADDBASE
    fprintf(stderr, "addBase() %c\n", ch);
#endif

    if (_style == 0)
      return(addBaseContiguous(cf, cr));

    if (_style == 1)
      return(addBaseCompressed(cf, cr));

    if (_style == 2)
      return(addBaseSpaced(cf, cr));

    if (_style == 3)
      return(addBaseCompressedSpaced(cf, cr));

    fprintf(stderr, "kMerBuilder::addBase()--  Invalid mer type %d.\n", _style);
    exit(1);

    return(false);
  }

  void    mask(void) {
    _fMer->mask(true);
    _rMer->mask(false);
  };

  kMer const   &theFMer(void) { return(*_fMer); };
  kMer const   &theRMer(void) { return(*_rMer); };
  kMer const   &theCMer(void) { return((theFMer() < theRMer()) ? theFMer() : theRMer()); };

  uint32        merSize(void)      { return(_merSize); };
  uint32        templateSpan(void) { return(_templateSpan); };

  uint32        baseSpan(uint32 b) {
    return(_compressionLength[(_compressionIndex + 1 + b) % _merSize]);;
  };

private:

  //  Style of builder we are
  uint32   _style;

  //  Amount of the mer that has valid sequence.  Sigh.  I really needed a signed value here --
  //  where negative values mean that we first have to get to the end of the template that was
  //  invalid, then we need to build a new mer.
  //
  //  And, yes, just simply making it signed leads to all sortes of compiler warnings about
  //  comparing signed and unsigned.  And I've been here before, and those warnings just propate
  //  endlessly.  Just go away, Mr. Smartypants.
  //
  //  Details: when building spaced seeds, if we hit an N in the middle of the template, we need to
  //  invalidate the mer, but not start building a new mer until we exhaust the current template.
  //  The example is template=1101.  Suppose we hit an N at the second 1.  We set the merSizeValid
  //  to 0, and proceed.  When we push on the base for the last 1 in the template, we'd increment
  //  the merSizeValid.  The first two 1's in the template would now create a mer big enough to be
  //  valid, and we'd return it -- but now the template we're using is 0111.
  //
  //  _merSizeValid is offset by _merSize (e.g., the true valid size is _merSizeValid - _merSize).
  //  _merSizeValidIs is the size _merSizeValid needs to be in order for it to be valid.
  //  Similarily, _merSizeValidZero is the value of zero (currently this is equal to _merSize).
  //
  uint32   _merSize;              //  Desired number of bases in the mer
  uint32  *_merSizeValid;         //  Actual number of bases in the mer
  uint32   _merSizeValidZero;     //  Definition of 'zero' bases in the mer
  uint32   _merSizeValidIs;       //  Definition of 'full' bases in the mer

  //  An array of mers, we allocate all mers in one block
  kMer    *_merStorage;

  //  Pointer to the currently active mer
  kMer    *_fMer;
  kMer    *_rMer;

  //  For compression
  uint32   _compression;
  uint32   _compressionIndex;        //  index into cL[] that is the last base in the mer
  uint32   _compressionFirstIndex;   //  index into cL[] that is the first base in a run
  uint32  *_compressionLength;       //  one per base
  uint32   _compressionCurrentLength;

  //  For templates
  uint32   _templateSpan;     //  # of 0's and 1's in the template
  uint32   _templateLength;   //  length of the pattern in the template
  char    *_template;         //  character string template
  uint32   _templatePos;      //  position we are building in the template
  uint32   _templateMer;      //  the mer we should output next
  uint32   _templateFirst;    //  if true, we're still building the initial mer
};

#endif  //  BIO_KMER_H
