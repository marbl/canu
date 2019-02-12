
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kmers.H"
#include "bits.H"

#include "files.H"


kmerCountStatistics::kmerCountStatistics() {
  _numUnique     = 0;
  _numDistinct   = 0;
  _numTotal      = 0;

  _histMax       = 32 * 1024 * 1024;      //  256 MB of histogram data.
  _hist          = new uint64 [_histMax];

  for (uint64 ii=0; ii<_histMax; ii++)
    _hist[ii] = 0;

  _histLen       = 0;
  _histVs        = NULL;
  _histOs        = NULL;
}


kmerCountStatistics::~kmerCountStatistics() {
  delete [] _hist;
  delete [] _histVs;
  delete [] _histOs;
}



void
kmerCountStatistics::clear(void) {
  _numUnique     = 0;
  _numDistinct   = 0;
  _numTotal      = 0;

  for (uint64 ii=0; ii<_histMax; ii++)
    _hist[ii] = 0;

  _histBig.clear();

  delete [] _histVs;
  delete [] _histOs;

  _histLen       = 0;
  _histVs        = NULL;
  _histOs        = NULL;
}



void
kmerCountStatistics::dump(stuffedBits *bits) {
  bits->setBinary(64, _numUnique);
  bits->setBinary(64, _numDistinct);
  bits->setBinary(64, _numTotal);

  //  Find the last used histogram value.

  uint32  histLast = _histMax;

  while (histLast-- > 0)
    if (_hist[histLast] > 0)
      break;

  histLast++;

  bits->setBinary(32, histLast);                //  Out of order relative to struct to keep
  bits->setBinary(32, _histBig.size());         //  the arrays below word-aligned.

  bits->setBinary(64, histLast, _hist);         //  Easy one, just dump an array.

  //  The histBig map is a little trickier.  The original format wrote two arrays,
  //  _hbigCount then _hbigNumber.  But those were always empty and never used.
  //  We're therefore free to change the order.

  for (map<uint64,uint64>::iterator it=_histBig.begin(); it != _histBig.end(); it++) {
    bits->setBinary(64, it->first);     //  Value
    bits->setBinary(64, it->second);    //  Number of occurrences
  }
}


void
kmerCountStatistics::dump(FILE        *outFile) {
  stuffedBits  *bits = new stuffedBits;

  dump(bits);

  bits->dumpToFile(outFile);

  delete bits;
}


void
kmerCountStatistics::load(stuffedBits *bits) {
  uint32  histLast;
  uint32  histBigLen;

  _numUnique   = bits->getBinary(64);
  _numDistinct = bits->getBinary(64);
  _numTotal    = bits->getBinary(64);

  histLast     = bits->getBinary(32);
  histBigLen   = bits->getBinary(32);

  delete [] _hist;

  _hist        = new uint64 [histLast + 1];
  _hist        = bits->getBinary(64, histLast, _hist);

  //  Count how many non-zero histogram values there are.

  _histLen = histBigLen;

  for (uint32 ii=0; ii<histLast; ii++)
    if (_hist[ii] > 0)
      _histLen++;

  //  Allocate space for them.

  _histVs = new uint32 [_histLen];
  _histOs = new uint64 [_histLen];

  //  Set the ones we've loaded already.

  _histLen = 0;

  for (uint32 ii=0; ii<histLast; ii++)
    if (_hist[ii] > 0) {
      _histVs[_histLen] = ii;
      _histOs[_histLen] = _hist[ii];
      _histLen++;
    }

  //  Load the rest from disk.

  for (uint32 ii=0; ii<histBigLen; ii++) {
    _histVs[_histLen] = bits->getBinary(64);
    _histOs[_histLen] = bits->getBinary(64);
    _histLen++;
  }

  //  Delete _hist to indicate we cannot accept new values.

  delete [] _hist;
  _hist = NULL;
}


void
kmerCountStatistics::load(FILE        *inFile) {
  stuffedBits  *bits = new stuffedBits;

  bits->loadFromFile(inFile);

  load(bits);

  delete bits;
}

