
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
 *    kmer/libbio/kmer.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-MAR-30 to 2014-APR-11
 *      are Copyright 2011,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kMer.H"

kMerBuilder::kMerBuilder(uint32 ms, uint32 cm, char *tm) {
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
    _compressionLength        = 0L;
    _compressionCurrentLength = 0;

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

    uint32 i=0, t=0;
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


  if (_merSize > KMER_WORDS * 32)
    fprintf(stderr, "kMer size too large; increase KMER_WORDS in libbio/kmer.H\n"), exit(1);


  _compressionLength = new uint32 [_merSize];

  for (uint32 z=0; z<_merSize; z++)
    _compressionLength[z] = (cm) ? 0 : 1;



  if (tm) {
    _merStorage   = new kMer   [_templateLength * 2];
    _merSizeValid = new uint32 [_templateLength];

    for (uint32 i=0; i<2*_templateLength; i++) {
      _merStorage[i].setMerSize(_merSize);
      _merStorage[i].setMerSpan(_templateSpan);
    }

    //  VERY IMPORTANT!  Offset the valid length to adjust for the
    //  template that every mer except the first is starting in the
    //  middle of.
    //
    for (uint32 i=0; i<_templateLength; i++)
      _merSizeValid[i] = _merSize - i;
  } else {
    _merStorage   = new kMer   [2];
    _merSizeValid = new uint32 [1];

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
}



kMerBuilder::~kMerBuilder() {
  delete [] _merSizeValid;
  delete [] _merStorage;
  delete [] _compressionLength;
  delete [] _template;
}



void
kMerBuilder::clear(bool clearMer) {

  //  Contiguous mers
  _merSizeValid[0] = _merSizeValidZero;

  //  Compressed mers
  if (_compression) {
    _compressionIndex         = 0;
    _compressionFirstIndex    = 0;
    _compressionCurrentLength = 0;

    for (uint32 z=0; z<_merSize; z++)
      _compressionLength[z] = 0;

    _merStorage[0].setMerSpan(0);
    _merStorage[1].setMerSpan(0);
  }

  //  Spaced mers
  if (_template) {
    for (uint32 i=0; i<2*_templateLength; i++)
      _merStorage[i].clear();

    for (uint32 i=0; i<_templateLength; i++)
      _merSizeValid[i] = _merSize - i;

    _templatePos     = 0;
    _templateMer     = 0;
    _templateFirst   = 1;
  }

  if (clearMer) {
    _fMer->clear();
    _rMer->clear();
  }
}




//
//  The addBase methods add a single base (cf - forward, cr - complemented) to
//  the mer.  The return true if another base is needed to finish the mer, and
//  false if the mer is complete.
//





bool
kMerBuilder::addBaseContiguous(uint64 cf, uint64 cr) {

  //  Not a valid base, reset the mer to empty, and request more bases
  //  (this is a slightly optimized version of clear()).
  if (cf & (unsigned char)0xfc) {
    clear(false);
    //_merSizeValid[0] = _merSizeValidZero;
    return(true);
  }

  //  Add the base to both mers.
  *_fMer += cf;
  *_rMer -= cr;

  //  If there aren't enough bases, request another one.
  if (_merSizeValid[0] + 1 < _merSizeValidIs) {
    _merSizeValid[0]++;
    return(true);
  }

  return(false);  //  Good!  Don't need another letter.
}






bool
kMerBuilder::addBaseCompressed(uint64 cf, uint64 cr) {

  //  Not a valid base, reset the mer to empty, and request more bases.
  //
  if (cf & (unsigned char)0xfc) {
    clear();
    return(true);
  }

  uint64  lb = theFMer().endOfMer(2);   //  Last base in the mer
  uint32  ms = theFMer().getMerSpan();  //  Span BEFORE adding the mer

  if (_merSizeValid[0] <= _merSizeValidZero)
    lb = 9;  //  No valid last base (should probably be ~uint64ZERO, but that screws up diagnostic output)

#ifdef DEBUGCOMP
  fprintf(stderr, "kMerBuilder::addBaseCompressed()--  lb="uint64FMT" cf="uint64FMT" ms="F_U32" ccl="F_U32" lvl="F_U32"\n",
          lb, cf, ms, _compressionCurrentLength, _compression);
#endif

  //  Always add one to the current length.  When we started, it
  //  was 0.  This represents the length AFTER adding the base.
  //
  _compressionCurrentLength++;

  //  If the lastbase is the same as the one we want to add (and
  //  there IS a last base), and we've seen too many of these,
  //  remember we've seen another letter in the run, and don't add
  //  it.  Request another letter.
  //
  if ((lb == cf) &&                                  //  last is the same as this
      (_compressionCurrentLength > _compression)) {  //  run is already too big
    _compressionLength[_compressionIndex]++;

    _fMer->setMerSpan(ms + 1);
    _rMer->setMerSpan(ms + 1);

#ifdef DEBUGCOMP
    fprintf(stderr, "kMerBuilder::addBaseCompressed()--  COMPRESSED currentIdx=%u first=%u",
            _compressionIndex, _compressionFirstIndex);
    for (uint32 x=0, y=_compressionFirstIndex; x<_merSize; x++) {
      fprintf(stderr, " %u(%d)", _compressionLength[y], y);
      y = (y + 1) % _merSize;
    }
    fprintf(stderr, "\n");
#endif
    return(true);
  }

  //  Else, it's a new run (a different letter) or our run isn't
  //  big enough to compress and we need to add the duplicate
  //  letter.

  *_fMer += cf;
  *_rMer -= cr;

  //  If this is a new letter, propagate the current length to the first letter in this run.  That
  //  way, when that letter is popped off the mer, we automagically update our span to include only
  //  as many letters as are here.
  //
  //                     01234567890
  //
  //  E.g.  For sequence TATTTTTTAGT (that's 6 T's) with a mersize of 3 and compression 2, we'd have
  //  mers with position:
  //
  //                                                            TATTTTTTAGT
  //  #1 TAT position 0 (with lengths 1, 1, 1) uncompressed mer TAT
  //  #2 ATT position 1 (with lengths 1, 1, 1)                   ATT
  //  #3 TTA position 6 (with lengths 5, 1, 1)                    TTTTTTA
  //  #4 TAG position 7                                                TAG
  //  #5 AGT position 8                                                 AGT
  //
  //  In #2, because the length so far (1) is not >= the compression (2) we add a new base and
  //  return.
  //
  //  In #3, the current length is >= the compression, so we keep stuffing on T's and incrementing
  //  the last length, stopping when we get the A.  We now propagate the current length to the first
  //  letter in the run.  Special case, if the first letter in the run is the first letter in the
  //  mer, we need to immediately update the span.

#ifdef DEBUGCOMP
  fprintf(stderr, "kMerBuilder::addBaseCompressed()--  ADDNEWBASE currentIdx=%u first=%u",
          _compressionIndex, _compressionFirstIndex);
  for (uint32 x=0, y=_compressionFirstIndex; x<_merSize; x++) {
    fprintf(stderr, " %u(%d)", _compressionLength[y], y);
    y = (y + 1) % _merSize;
  }
  fprintf(stderr, "\n");
#endif

  //  If we added a new letter, transfer the run-length count to the first letter in the previous
  //  run.  In the above example, when we built the run, the lengths are (1, 1, 5).  That is, all
  //  compression occurred on the last letter.  When we shift off that first letter, we want to
  //  remove as much of the run as possible.

  if (lb != cf) {
    if (_compressionFirstIndex != _compressionIndex) {
      _compressionLength[_compressionFirstIndex] += _compressionLength[_compressionIndex] - 1;
      _compressionLength[_compressionIndex]       = 1;
    }
    _compressionFirstIndex    = (_compressionIndex + 1) % _merSize;
    _compressionCurrentLength = 1;
  }

  _compressionIndex = (_compressionIndex + 1) % _merSize;
  ms -= _compressionLength[_compressionIndex];  //  subtract the count for the letter we just shifted out

#ifdef DEBUGCOMP
  fprintf(stderr, "kMerBuilder::addBaseCompressed()--  ADDNEWBASE shifted out at idx="F_U32" with "F_U32" positions; final span "F_U32"\n",
          _compressionIndex,
          _compressionLength[_compressionIndex],
          ms + 1);
#endif

  _compressionLength[_compressionIndex] = 1;    //  one letter at this position

  _fMer->setMerSpan(ms + 1);
  _rMer->setMerSpan(ms + 1);

  //  If there aren't enough bases, request another one.
  if (_merSizeValid[0] + 1 < _merSizeValidIs) {
    _merSizeValid[0]++;
    return(true);
  }

  return(false);  //  Good!  Don't need another letter.
}






bool
kMerBuilder::addBaseSpaced(uint64 cf, uint64 cr) {
#ifdef DEBUGSPACE
  fprintf(stderr, "add %c templatePos=%u templateMer=%u\n", ch, _templatePos, _templateMer);
#endif

  //  We always advance the templatePos, unfortunately, we need to
  //  use the current value throughout this function.  If there
  //  was a single return point, we could advance immediately
  //  before returning.
  //
  uint32 tp = _templatePos;
  _templatePos = (_templatePos + 1) % _templateLength;

  //  If we get an invalid letter, set all mers that would have
  //  had a letter added to be broken.
  //
  if (cf & (unsigned char)0xfc) {

    for (uint32 m=0; m<_templateLength; m++) {
      uint32 tppos = (tp + _templateLength - m) % _templateLength;

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
  for (uint32 m=0; m<_templateLength; m++) {
    uint32  tppos = (tp + _templateLength - m) % _templateLength;

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






bool
kMerBuilder::addBaseCompressedSpaced(uint64 UNUSED(cf), uint64 UNUSED(cr)) {
  fprintf(stderr, "kMerBuilder::addBaseCompressedSpace()--  Compressed and spaced mers not supported.\n");
  exit(1);
}

