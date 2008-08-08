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


#undef DEBUG
#undef DEBUGCOMP
#undef DEBUGSPACE


class kMerBuilder {
public:
  kMerBuilder(u32bit ms=0, u32bit cm=0, char *tm=0L) {
    _style            = 0;

    _merSize          = 0;
    _merSizeValid     = 0L;
    _merSizeValidIs   = 0;
    _merSizeValidZero = 0;

    _merStorage       = 0L;
    _fMer             = 0L;
    _rMer             = 0L;

   _compression               = 0;
   _compressionIndex          = 0;
   _compressionFirstIndex     = 0;
   _compressionLength         = 0L;
   _compressionCurrentLength  = 0;

   _templateSpan    = 0;
   _templateLength  = 0;
   _template        = 0L;
   _templatePos     = 0;
   _templateMer     = 0;
   _templateFirst   = 0;

   if (ms) {
     _style            = 0;
     _merSize          = ms;
     _merSizeValidIs   = _merSize + _merSize;
     _merSizeValidZero = _merSize;
   }

   if (cm) {
     _style = 1;

     _merSize          = ms;
     _merSizeValidIs   = _merSize + _merSize;
     _merSizeValidZero = _merSize;

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
     _style           = 2;
     _merSize         = 0;
     _templateSpan    = strlen(tm);
     _templateLength  = 0;
     _template        = new char [_templateSpan + 1];
     _templatePos     = 0;
     _templateMer     = 0;
     _templateFirst   = 1;

     //  Templates cannot begin or end in zero -- they MUST begin/end
     //  with a letter.  We silently fix these problems.  Unless there
     //  are no 1's in the string, then we bail.

     u32bit i=0, t=0;
     while ((i < _templateSpan) && (tm[i] == '0'))
       i++;

     if (i == _templateSpan) {
       fprintf(stderr, "invalid kMerBuilder template '%s' -- its empty!\n", tm);
       exit(1);
     }

     while (i < _templateSpan) {
       _template[t] = 0;
       
       if (tm[i] == '1') {
         _template[t] = 1;
         _merSize++;
       }

       i++;
       t++;
     }

     while (_template[--t] == 0)
       ;

     _templateSpan = t + 1;

#ifdef DEBUGSPACE
     for (i=0; i<_templateSpan; i++)
       fprintf(stderr, "%d", _template[i]);
     fprintf(stderr, " -- %d\n", _templateSpan);
#endif

     //  Look for patterns in the template, set templateLength to be the
     //  size of the pattern.

     _templateLength = _templateSpan;

     //  Finally, we can set what valid and zero mersizes are.
     _merSizeValidIs   = _templateLength + _merSize;
     _merSizeValidZero = _templateLength;
   }

   if (cm && tm) {
     _style = 3;
     assert(0);
   }

   if (tm) {
     _merStorage   = new kMer   [_templateLength * 2];
     _merSizeValid = new u32bit [_templateLength];

     for (u32bit i=0; i<2*_templateLength; i++) {
       _merStorage[i].setMerSize(_merSize);
       _merStorage[i].setMerSpan(_templateSpan);
     }

     //  VERY IMPORTANT!  Offset the valid length to adjust for the
     //  template that every mer except the first is starting in the
     //  middle of.
     //
     for (u32bit i=0; i<_templateLength; i++)
       _merSizeValid[i] = _merSize - i;
   } else {
     _merStorage   = new kMer   [2];
     _merSizeValid = new u32bit [1];

     _merStorage[0].setMerSize(_merSize);
     _merStorage[1].setMerSize(_merSize);

     _merSizeValid[0] = _merSizeValidZero;

     if (cm) {
       _merStorage[0].setMerSpan(0);
       _merStorage[1].setMerSpan(0);
     }
   }

   _fMer = _merStorage + 0;
   _rMer = _merStorage + 1;
  };


  ~kMerBuilder() {
    delete [] _merStorage;
    delete [] _merSizeValid;
  };


  //  Clear all mer data
  void    clear(void) {

    if (_template) {
      for (u32bit i=0; i<2*_templateLength; i++)
        _merStorage[i].clear();
     for (u32bit i=0; i<_templateLength; i++)
       _merSizeValid[i] = _merSize - i;
     _templatePos     = 0;
     _templateMer     = 0;
     _templateFirst   = 1;
    }

    if (_compression) {
      _compressionIndex         = 0;
      _compressionFirstIndex    = 0;
      _compressionCurrentLength = 0;

      for (u32bit z=0; z<_merSize; z++)
        _compressionLength[z] = 0;

      _merStorage[0].setMerSpan(0);
      _merStorage[1].setMerSpan(0);
    }

    _fMer->clear();
    _rMer->clear();
    _merSizeValid[0] = _merSizeValidZero;
  };

  //  Returns true if we need another base to finish the mer.  This
  //  only occurs for compressed mers, if we are in a homopolymer run.
  //
  bool    addBase(char ch) {
    u64bit  cf = letterToBits[ch];
    u64bit  cr = letterToBits[complementSymbol[ch]];

#ifdef DEBUG
    fprintf(stderr, "addBase() %c\n", ch);
#endif

    ////////////////////////////////////////
    //
    //  Normal contiguous mers
    //
    if (_style == 0) {

      if (cf & (unsigned char)0xfc) {
        _merSizeValid[0] = _merSizeValidZero;
        return(true);
      }

      *_fMer += cf;
      *_rMer -= cr;

      if (_merSizeValid[0] + 1 < _merSizeValidIs) {
        _merSizeValid[0]++;
        return(true);
      }
      return(false);  //  Good!  Don't need another letter.
    }

    ////////////////////////////////////////
    //
    //  Compressed mers
    //
    if (_style == 1) {

      //  If it's an invalid letter, get the heck outta here, reset stuff
      //
      //  Unfortunately, we still need to set all the lengths to zero.
      //  The next time after an invalid letter we'd be subtracting
      //  out a bogus value.
      //
      if (cf & (unsigned char)0xfc) {
        _compressionIndex         = 0;  //_merSize - 1;  //  Next time, this will be incremented to 0
        _compressionFirstIndex    = 0;
        _compressionCurrentLength = 0;

        for (u32bit z=0; z<_merSize; z++)
          _compressionLength[z] = 0;

        _fMer->setMerSpan(0);
        _rMer->setMerSpan(0);

        _merSizeValid[0] = _merSizeValidZero;

#ifdef DEBUGCOMP
        fprintf(stderr, "invalid letter, keep going.\n");
#endif
        return(true);
      }

      //  Else, valid letter.

      u64bit  lb = theFMer().endOfMer(2);
      u32bit  ms = theFMer().getMerSpan();

#ifdef DEBUGCOMP
      fprintf(stderr, "lb="u64bitFMT" cf="u64bitFMT" ms="u32bitFMT" ccl="u32bitFMT" compression="u32bitFMT"\n",
              lb, cf, ms, _compressionCurrentLength, _compression);
#endif

      //  Always add one to the current length.  When we started, it
      //  was 0, and this represents the length after adding the base.
      //
      _compressionCurrentLength++;

      //  If the lastbase is the same as the one we want to add (and
      //  there IS a last base), and we've seen too many of these,
      //  remember we've seen another letter in the run, and don't add
      //  it.  Request another letter.
      //
      if ((_merSizeValidZero < _merSizeValid[0]) &&      //  there is a last base
          (lb == cf) &&                                  //  last is the same as this
          (_compressionCurrentLength > _compression)) {  //  run is already too big
        _compressionLength[_compressionIndex]++;

        _fMer->setMerSpan(ms + 1);
        _rMer->setMerSpan(ms + 1);

#ifdef DEBUGCOMP
        {
          for (int x=0, y=_compressionFirstIndex; x<_merSize; x++) {
            fprintf(stderr, "%u(%d) ", _compressionLength[y], y);
            y = (y + 1) % _merSize;
          }
          fprintf(stderr, " currentIdx=%u first=%u\n", _compressionIndex, _compressionFirstIndex);
        }
        fprintf(stderr, "need more!\n");
#endif
        return(true);
      } else {
#ifdef DEBUGCOMP
        fprintf(stderr, "msvis="u64bitFMT" <= msv="u64bitFMT"\n",
                _merSizeValidIs, _merSizeValid[0]);
#endif
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

#ifdef DEBUGCOMP
      {
        for (int x=0, y=_compressionFirstIndex; x<_merSize; x++) {
          fprintf(stderr, "%u(%d) ", _compressionLength[y], y);
          y = (y + 1) % _merSize;
        }
        fprintf(stderr, " currentIdx=%u first=%u MID\n", _compressionIndex, _compressionFirstIndex);
      }
#endif

      //  If it's a new letter, move the count to the first thing in
      //  the run.
      //
      if (lb != cf) {
        if (_compressionFirstIndex != _compressionIndex) {
          _compressionLength[_compressionFirstIndex] += _compressionLength[_compressionIndex] - 1;
          _compressionLength[_compressionIndex]       = 1;
        }
        _compressionFirstIndex    = (_compressionIndex + 1) % _merSize;
        _compressionCurrentLength = 1;
      } else {
      }

      _compressionIndex = (_compressionIndex + 1) % _merSize;
      ms -= _compressionLength[_compressionIndex];  //  subtract the count for the letter we just shifted out
      _compressionLength[_compressionIndex] = 1;    //  one letter at this position

      _fMer->setMerSpan(ms + 1);
      _rMer->setMerSpan(ms + 1);

      if (_merSizeValid[0] + 1 < _merSizeValidIs) {
        _merSizeValid[0]++;
#ifdef DEBUGCOMP
        {
          for (int x=0, y=_compressionFirstIndex; x<_merSize; x++) {
            fprintf(stderr, "%d ", _compressionLength[y]);
            y = (y + 1) % _merSize;
          }
          fprintf(stderr, " currentIdx=%u first=%u\n", _compressionIndex, _compressionFirstIndex);
        }
        fprintf(stderr, "still too short, need more!\n");
#endif
        return(true);
      }

      //  Good!  Don't need another letter.
#ifdef DEBUGCOMP
        {
          for (int x=0, y=_compressionFirstIndex; x<_merSize; x++) {
            fprintf(stderr, "%d ", _compressionLength[y]);
            y = (y + 1) % _merSize;
          }
          fprintf(stderr, " currentIdx=%u first=%u\n", _compressionIndex, _compressionFirstIndex);
        }
      fprintf(stderr, "span=%u %u\n", ms+1, _fMer->getMerSpan());
#endif
      return(false);
    }

    ////////////////////////////////////////
    //
    //  Spaced mers
    //
    if (_style == 2) {

#ifdef DEBUGSPACE
      fprintf(stderr, "add %c templatePos=%u templateMer=%u\n", ch, _templatePos, _templateMer);
#endif

      //  We always advance the templatePos, unfortunately, we need to
      //  use the current value throughout this function.  If there
      //  was a single return point, we could advance immediately
      //  before returning.
      //
      u32bit tp = _templatePos;
      _templatePos = (_templatePos + 1) % _templateLength;

      //  If we get an invalid letter, set all mers that would have
      //  had a letter added to be broken.
      //
      if (cf & (unsigned char)0xfc) {

        for (u32bit m=0; m<_templateLength; m++) {
          u32bit tppos = (tp + _templateLength - m) % _templateLength;

          if (_template[tppos] == 1) {

            //  Reset to 'zero', but make it skip over any remaining
            //  positions in the current template.
            //
            _merSizeValid[m]  = _merSizeValidZero + tppos - _templateLength + 1;

#ifdef DEBUGSPACE
            fprintf(stderr, "-- invalid letter, reset mer %u to valid %u (mersizevalidzero=%u ttpos=%u templatelength=%u)\n",
                    m, _merSizeValid[m], _merSizeValidZero, tppos, _templateLength);
#endif
          }
        }

        if (_templateFirst == 0)
          _templateMer = (_templateMer + 1) % _templateLength;

        return(true);
      }

      //  We have a valid letter, and add it to all the mers that the
      //  template allows.
      //
      for (u32bit m=0; m<_templateLength; m++) {
        u32bit  tppos = (tp + _templateLength - m) % _templateLength;

        if (_template[tppos] == 1) {
          _merStorage[2*m+0] += cf;
          _merStorage[2*m+1] -= cr;

          if (_merSizeValid[m] < _merSizeValidIs)
            _merSizeValid[m]++;

#ifdef DEBUGSPACE
          fprintf(stderr, "push %c onto %d (at template %u)  length = %u  %s\n",
                  ch, m, (tp + _templateLength - m) % _templateLength,
                  _merSizeValid[m],
                  (_merSizeValid[m] >= _merSizeValidIs) ? "complete" : "");
#endif
        } else if (_merSizeValid[m] <= _merSizeValidZero) {

          //  The template doesn't want us to add a letter to the mer,
          //  but we're adjusting for an aborted template, and we're
          //  counting template positions (not just non-zero template
          //  positions) when adjusting.
          //
          _merSizeValid[m]++;
        }
      }

      //  If the current mer isn't long enough, we move to the next mer,
      //  and request another letter.
      //
      if (_merSizeValid[_templateMer] < _merSizeValidIs) {
        if (_templateFirst == 0)
          _templateMer = (_templateMer + 1) % _templateLength;
#ifdef DEBUGSPACE
        fprintf(stderr, "-- too short -- need more templateMer=%u templateFirst=%u\n", _templateMer, _templateFirst);
#endif
        return(true);
      }

      //  On startup, _templateMer is always 0 (the first mer) until
      //  it is long enough to be a valid mer.  Then, we clear
      //  _templateFirst so that we can start advancing through mers.

      //  Update the f and r pointers to the correct mers, advance our
      //  template to the next, and terminate.
      //
      _fMer = _merStorage + 2 * _templateMer + 0;
      _rMer = _merStorage + 2 * _templateMer + 1;

#ifdef DEBUGSPACE
      fprintf(stderr, "-- valid!  (templateMer = %u)\n", _templateMer);
#endif

      _templateFirst = 0;
      _templateMer   = (_templateMer + 1) % _templateLength;

      return(false);  //  Good!  Don't need another letter.
    }

    ////////////////////////////////////////
    //
    //  Spaced mers, using compression
    //
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

  //  Style of builder we are
  u32bit   _style;

  //  Amount of the mer that has valid sequence.  Sigh.  I really
  //  needed a signed value here -- where negative values mean that we
  //  first have to get to the end of the template that was invalid,
  //  then we need to build a new mer.
  //
  //  And, yes, just simply making it signed leads to all sortes of
  //  compiler warnings about comparing signed and unsigned.  And I've
  //  been here before, and those warnings just propate endlessly.
  //  Just go away, Mr. Smartypants.
  //
  //  Details: when building spaced seeds, if we hit an N in the
  //  middle of the template, we need to invalidate the mer, but not
  //  start building a new mer until we exhaust the current template.
  //  The example is template=1101.  Suppose we hit an N at the second
  //  1.  We set the merSizeValid to 0, and proceed.  When we push on
  //  the base for the last 1 in the template, we'd increment the
  //  merSizeValid.  The first two 1's in the template would now
  //  create a mer big enough to be valid, and we'd return it -- but
  //  now the template we're using is 0111.
  //
  //  _merSizeValid is offset by _merSize (e.g., the true valid size
  //  is _merSizeValid - _merSize).  _merSizeValidIs is the valid
  //  _merSizeValid needs to be in order for it to be valid.
  //  Similarily, _merSizeValidZero is the value of zero (currently
  //  this is equal to _merSize).
  //
  u32bit   _merSize;              //  Number of bases in the mer
  u32bit  *_merSizeValid;
  u32bit   _merSizeValidZero;
  u32bit   _merSizeValidIs;

  //  An array of mers, we allocate all mers in one block
  kMer    *_merStorage;

  //  Pointer to the currently active mer
  kMer    *_fMer;
  kMer    *_rMer;

  //  For compression
  u32bit   _compression;
  u32bit   _compressionIndex;   //  index into cL[] that is the last base in the mer
  u32bit   _compressionFirstIndex;
  u32bit  *_compressionLength;  //  one per base
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
