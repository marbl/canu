
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

  //  This only writes the latest version.

  bits->setBinary(64, _numUnique);
  bits->setBinary(64, _numDistinct);
  bits->setBinary(64, _numTotal);

  //  Find the maximum value, and count how many values we have
  //  in the histogram.

  uint64   numValues = _histBig.size();

  for (uint32 ii=0; ii<_histMax; ii++)
    if (_hist[ii] > 0)
      numValues++;

  //

  bits->setBinary(64, numValues);

  //  Now the data!

  for (uint32 ii=0; ii<_histMax; ii++) {
    if (_hist[ii] > 0) {
      bits->setBinary(64,       ii);     //  Value
      bits->setBinary(64, _hist[ii]);    //  Number of occurrences
    }
  }

  for (map<uint64,uint64>::iterator it=_histBig.begin(); it != _histBig.end(); it++) {
    bits->setBinary(64, it->first);      //  Value
    bits->setBinary(64, it->second);     //  Number of occurrences
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
kmerCountStatistics::load_v01(stuffedBits *bits) {
  uint32  histLast;
  uint32  hbigLen;

  _numUnique   = bits->getBinary(64);
  _numDistinct = bits->getBinary(64);
  _numTotal    = bits->getBinary(64);

  histLast     = bits->getBinary(32);
  hbigLen      = bits->getBinary(32);

  //fprintf(stderr, "kmerCountStatistics::load_v01()-- %lu %lu %lu %u %u\n",
  //        _numUnique, _numDistinct, _numTotal, histLast, hbigLen);

  //  Load the histogram values.

  uint64  *hist = new uint64 [histLast + 1];

  hist         = bits->getBinary(64, histLast, hist);

  //  (over) allocate space for the histogram list.

  _histVs = new uint64 [histLast + hbigLen + 1];
  _histOs = new uint64 [histLast + hbigLen + 1];

  //  Convert the loaded hist[] into _histVs and _histOs.

  _histLen = 0;

  for (uint32 ii=0; ii<histLast; ii++) {
    if (_hist[ii] > 0) {
      _histVs[_histLen] = ii;
      _histOs[_histLen] = hist[ii];
      _histLen++;
    }
  }

  delete [] hist;

#if 0
  //  If hbigLen isn't zero, we have the intermediate format, that lived for
  //  about a day, that stores large values too.

  if (hbigLen > 0) {
    for (uint64 ii=0; ii<hbigLen; ii++) {
      _histVs[_histLen] = bits->getBinary(64);
      _histOs[_histLen] = bits->getBinary(64);
      _histLen++;
    }
  }
#endif

  //  Delete _hist to indicate we cannot accept new values.

  delete [] _hist;
  _hist = NULL;
}



void
kmerCountStatistics::load_v03(stuffedBits *bits) {

  _numUnique   = bits->getBinary(64);
  _numDistinct = bits->getBinary(64);
  _numTotal    = bits->getBinary(64);
  _histLen     = bits->getBinary(64);

  //fprintf(stderr, "kmerCountStatistics::load_v03()-- %lu %lu %lu %lu\n",
  //        _numUnique, _numDistinct, _numTotal, _histLen);

  //  Allocate space.

  _histVs = new uint64 [_histLen];
  _histOs = new uint64 [_histLen];

  //  Load the values into our list.

  for (uint64 ii=0; ii<_histLen; ii++) {
    _histVs[ii] = bits->getBinary(64);
    _histOs[ii] = bits->getBinary(64);
  }

  //  Delete _hist to indicate we cannot accept new values.

  delete [] _hist;
  _hist = NULL;
}



void
kmerCountStatistics::load(stuffedBits *bits,
                          uint32       version) {

  switch (version) {
    case 1:
    case 2:
      load_v01(bits);
      break;
    case 3:
      load_v03(bits);
      break;
    default:
      fprintf(stderr, "kmerCountStatistics::load()-- Unknown version %u\n", version), exit(1);
      break;
  }
}



void
kmerCountStatistics::load(FILE        *inFile,
                          uint32       version) {
  stuffedBits  *bits = new stuffedBits;

  bits->loadFromFile(inFile);

  load(bits, version);

  delete bits;
}

