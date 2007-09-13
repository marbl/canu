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

   _compression    = cm;

   _templateSize   = 0;
   _templateSpan   = 0;
   _template       = 0L;

   if (ms) {
     _style        = 0;
     _merSize      = ms;
   }

   if (cm) {
     _style        = 1;
     _compression  = cm;
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
     }
     _fMer = _fmers[_currentMer];
     _rMer = _rmers[_currentMer];
   } else {
     _merStorage = new kMer [2];

     _fMer = _merStorage + 0;
     _rMer = _merStorage + 1;

     _fMer->setMerSize(_merSize);
     _rMer->setMerSize(_merSize);
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
    bool    rr = false;

    switch (_style) {
      case 0:
        //  Normal contiguous mers
        *_fMer += cf;
        *_rMer -= cr;

        if (_merSizeValid < _merSize)
          _merSizeValid++;
        if (cf & (unsigned char)0xfc)
          _merSizeValid = 0;
        if (_merSizeValid < _merSize)
          rr = true;
        break;
      case 1:
        //  Compressed mers
        assert(0);
        break;
      case 2:
        //  Spaced mers
        assert(0);
        break;
      case 3:
        //  Spaced mers, using compression
        assert(0);
        break;
      default:
        assert(0);
        break;
    }

    return(rr);
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

  u32bit   _compression;

  u32bit   _templateSize;  //  # of 0's and 1's in the template
  u32bit   _templateSpan;  //  If the template is repetitive, this is the length of the pattern.
  char    *_template;
};

#endif  //  BIO_KMER_H
