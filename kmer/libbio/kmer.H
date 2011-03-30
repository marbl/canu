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

#include "util.h"
#include "util++.H"
#include "bio.h"
#include "bio++.H"
#include "kmeriface.H"

#if KMER_WORDS == 1
#include "kmertiny.H"
typedef kMerTiny kMer;
#else
#include "kmerhuge.H"
typedef kMerHuge kMer;
#endif


#undef DEBUGADDBASE
#undef DEBUGCOMP
#undef DEBUGSPACE


class kMerBuilder {
public:
  kMerBuilder(u32bit ms=0, u32bit cm=0, char *tm=0L);
  ~kMerBuilder();

  //  Clear all mer data, reset state to as just after construction.
  void    clear(bool clearMer=true);

  //  Returns true if we need another base to finish the mer.  This
  //  only occurs for compressed mers, if we are in a homopolymer run.
  //
private:
  bool addBaseContiguous(u64bit cf, u64bit cr);
  bool addBaseCompressed(u64bit cf, u64bit cr);
  bool addBaseSpaced(u64bit cf, u64bit cr);
  bool addBaseCompressedSpaced(u64bit cf, u64bit cr);

public:
  bool    addBase(char ch) {
    u64bit  cf = letterToBits[ch];
    u64bit  cr = letterToBits[complementSymbol[ch]];

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

  u32bit        merSize(void)      { return(_merSize); };
  u32bit        templateSpan(void) { return(_templateSpan); };

  u32bit        baseSpan(u32bit b) {
    return(_compressionLength[(_compressionIndex + 1 + b) % _merSize]);;
  };

private:

  //  Style of builder we are
  u32bit   _style;

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
  u32bit   _merSize;              //  Desired number of bases in the mer
  u32bit  *_merSizeValid;         //  Actual number of bases in the mer
  u32bit   _merSizeValidZero;     //  Definition of 'zero' bases in the mer
  u32bit   _merSizeValidIs;       //  Definition of 'full' bases in the mer

  //  An array of mers, we allocate all mers in one block
  kMer    *_merStorage;

  //  Pointer to the currently active mer
  kMer    *_fMer;
  kMer    *_rMer;

  //  For compression
  u32bit   _compression;
  u32bit   _compressionIndex;        //  index into cL[] that is the last base in the mer
  u32bit   _compressionFirstIndex;   //  index into cL[] that is the first base in a run
  u32bit  *_compressionLength;       //  one per base
  u32bit   _compressionCurrentLength;

  //  For templates
  u32bit   _templateSpan;     //  # of 0's and 1's in the template
  u32bit   _templateLength;   //  length of the pattern in the template
  char    *_template;         //  character string template
  u32bit   _templatePos;      //  position we are building in the template
  u32bit   _templateMer;      //  the mer we should output next
  u32bit   _templateFirst;    //  if true, we're still building the initial mer
};

#endif  //  BIO_KMER_H
