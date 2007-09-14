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

#include "kmeriface.H"

#if KMER_WORDS == 1
#include "kmertiny.H"
typedef kMerTiny kMer;
#else
#include "kmerhuge.H"
typedef kMerHuge kMer;
#endif




class kMerBuilder {
public:
  kMerBuilder(u32bit ms=0, u32bit cm=0, char *tm=0L) {
    _fMer          = 0L;
    _rMer          = 0L;

    _currentMer    = 0;
    _fmers         = 0L;
    _rmers         = 0L;

    _merSizeValid  = 0;

    _style         = 0;

    _merSize       = ms;

   _compression               = cm;
   _compressionIndex          = 0;
   _compressionFirstIndex     = 0;
   _compressionLength         = 0L;
   _compressionCurrentLength  = 0;

   _templateSize   = 0;
   _templateSpan   = 0;
   _template       = 0L;

   if (ms) {
     _style = 0;
     _merSize = ms;
   }

   if (cm) {
     _style = 1;

     _compression              = cm;
     _compressionIndex         = 0;
     _compressionFirstIndex    = 0;
     _compressionLength        = new u32bit [_merSize];
     _compressionCurrentLength = 0;

     for (u32bit z=0; z<_merSize; z++)
       _compressionLength[z] = 0;

     assert(_compression < _merSize);
   }

   if (tm) {
     _style        = 2;
     _merSize      = 0;
     _templateSize = strlen(tm);
     _templateSpan = 0;

     _template     = new char [_templateSpan + 1];
     strcpy(_template, tm);

     for (u32bit i=0; i<_templateSize; i++)
       if (_template[i] == '1')
         _merSize++;

     //  Look for patterns in the template, set templateSpan to be the
     //  size of the pattern.

     _templateSpan = _templateSize;
   }

   if (cm && tm) {
     _style = 3;
   }

   if (tm) {
     _merStorage = new kMer [2 * _templateSpan];

     for (u32bit i=0; i<_templateSpan; i++) {
       _fmers = new kMer* [_templateSpan];
       _rmers = new kMer* [_templateSpan];

       _fmers[i] = _merStorage + 2 * i;
       _rmers[i] = _merStorage + 2 * i + 1;

       _fmers[i]->setMerSize(_merSize);
       _rmers[i]->setMerSize(_merSize);

       _fmers[i]->setMerSpan(_templateSize);
       _rmers[i]->setMerSpan(_templateSize);
     }
     _fMer = _fmers[_currentMer];
     _rMer = _rmers[_currentMer];
   } else {
     _merStorage = new kMer [2];

     _fMer = _merStorage + 0;
     _rMer = _merStorage + 1;

     _fMer->setMerSize(_merSize);
     _rMer->setMerSize(_merSize);

     if (cm) {
       _fMer->setMerSpan(0);
       _rMer->setMerSpan(0);
     }
   }
  };


  ~kMerBuilder() {
    delete [] _merStorage;
  };


  //  Clear all mer data
  void    clear(void) {
    if (_template) {
      for (u32bit i=0; i<_templateSpan; i++) {
        _fmers[i]->clear();
        _rmers[i]->clear();
      }
    } else {
      _fMer->clear();
      _rMer->clear();
    }
    _merSizeValid = 0;
  };

  //  Returns true if we need another base to finish the mer.  This
  //  only occurs for compressed mers, if we are in a homopolymer run.
  //
  bool    addBase(char ch) {
    u64bit  cf = validCompressedSymbol[ch];
    u64bit  cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

    //  Normal contiguous mers
    if (_style == 0) {

      if (cf & (unsigned char)0xfc) {
        _merSizeValid = 0;
        return(true);
      }

      *_fMer += cf;
      *_rMer -= cr;

      if (_merSizeValid + 1 < _merSize) {
        _merSizeValid++;
        return(true);
      }
      return(false);  //  Good!  Don't need another letter.
    }


    //  Compressed mers
    if (_style == 1) {

      //  If it's an invalid letter, get the heck outta here, reset stuff
      //
      //  Unfortunately, we still need to set all the lengths to zero.
      //  The next time after an invalid letter we'd be subtracting
      //  out a bogus value.
      //
      if (cf & (unsigned char)0xfc) {
        _compressionIndex         = _merSize - 1;  //  Next time, this will be incremented to 0
        _compressionFirstIndex    = 0;
        _compressionCurrentLength = 0;

        for (u32bit z=0; z<_merSize; z++)
          _compressionLength[z] = 0;

        _fMer->setMerSpan(0);
        _rMer->setMerSpan(0);

        _merSizeValid = 0;

        return(true);
      }

      //  Else, valid letter.

      u64bit  lb = theFMer().endOfMer(2);
      u32bit  ms = theFMer().getMerSpan();

      //  Always add one to the current length.  When we started, it
      //  was 0, and this represents the length after adding the base.
      //
      _compressionCurrentLength++;

      //  If the lastbase is the same as the one we want to add (and
      //  there IS a last base), and we've seen too many of these,
      //  remember we've seen another letter in the run, and don't add
      //  it.  Request another letter.
      //
      if ((_merSizeValid > 0) &&
          (lb == cf) &&
          (_compressionCurrentLength > _compression)) {
        _compressionLength[_compressionIndex]++;

        _fMer->setMerSpan(ms + 1);
        _rMer->setMerSpan(ms + 1);

        return(true);
      }

      //  Else, it's a new run (a different letter) or our run isn't
      //  big enough to compress and we need to add the duplicate
      //  letter.

      *_fMer += cf;
      *_rMer -= cr;

      //  If this is a new letter, propagate the current length to the
      //  first letter in this run.  That way, when that letter is
      //  popped off the mer, we automagically update our span to
      //  include only as many letters as are here.
      //
      //                     01234567890
      //  E.g.  For sequence TATTTTTTAGT (that's 6 T's) with a
      //  mersize of 3 and compression 2, we'd have mers with
      //  position:
      //
      //  #1 TAT position 0 (with lengths 1, 1, 1)
      //  #2 ATT position 1 (with lengths 1, 1, 1)
      //  #3 TTA position 6 (with lengths 5, 1, 1)
      //  #4 TAG position 7
      //  #5 AGT position 8
      //
      //  In #2, because the length so far (1) is not >= the
      //  compression (2) we add a new base and return.
      //
      //  In #3, the current length is >= the compression, so we keep
      //  stuffing on T's and incrementing the last length, stopping
      //  when we get the A.  We now propagate the current length to
      //  the first letter in the run.  Special case, if the first
      //  letter in the run is the first letter in the mer, we need to
      //  immediately update the span.

      //  If it's a new letter, move the count to the first thing in
      //  the run.
      //
      if (lb != cf) {
        _compressionLength[_compressionFirstIndex] += _compressionLength[_compressionIndex] - 1;
        _compressionLength[_compressionIndex]       = 1;
        _compressionFirstIndex    = (_compressionIndex + 1) % _merSize;
        _compressionCurrentLength = 1;
      } else {
      }

      _compressionIndex = (_compressionIndex + 1) % _merSize;
      ms -= _compressionLength[_compressionIndex];  //  subtract the count for the letter we just shifted out
      _compressionLength[_compressionIndex] = 1;    //  one letter at this position

      _fMer->setMerSpan(ms + 1);
      _rMer->setMerSpan(ms + 1);

      if (_merSizeValid + 1 < _merSize) {
        _merSizeValid++;
        return(true);
      }

      //fprintf(stderr, "span=%u %u\n", ms+1, _fMer->getMerSpan());

      return(false);  //  Good!  Don't need another letter.
    }

    //  Spaced mers
    if (_style == 2) {
      assert(0);
    }

    //  Spaced mers, using compression
    if (_style == 3) {
      assert(0);
    }

    assert(0);
    return(false);
  }

  void    mask(void) {
    _fMer->mask(true);
    _rMer->mask(false);
  };

  kMer const   &theFMer(void) { return(*_fMer); };
  kMer const   &theRMer(void) { return(*_rMer); };
  kMer const   &theCMer(void) { return((theFMer() < theRMer()) ? theFMer() : theRMer()); };

private:

  //  An array of mers, we allocate all mers in one block
  kMer    *_merStorage;

  //  Pointer to the currently active mer
  kMer    *_fMer;
  kMer    *_rMer;

  //  If spaced, we need to keep a set of mers
  u32bit   _currentMer;
  kMer   **_fmers;
  kMer   **_rmers;

  //  Amount of the mer that has valid sequence
  u32bit   _merSizeValid;
  //u32bit  *_merSizeValid;

  u32bit   _style;

  u32bit   _merSize;       //  Number of bases in the mer

  //  For compression
  u32bit   _compression;
  u32bit   _compressionIndex;   //  index into cL[] that is the last base in the mer
  u32bit   _compressionFirstIndex;
  u32bit  *_compressionLength;  //  one per base
  u32bit   _compressionCurrentLength;

  //  For templates
  u32bit   _templateSize;  //  # of 0's and 1's in the template
  u32bit   _templateSpan;  //  If the template is repetitive, this is the length of the pattern.
  char    *_template;
};

#endif  //  BIO_KMER_H
