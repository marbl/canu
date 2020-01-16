
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

#ifndef LIBKMER_HISTOGRAM_H
#define LIBKMER_HISTOGRAM_H

#ifndef LIBKMER_H
#error "include kmers.H, not this."
#endif

#include "kmers-tiny.H"

#include <map>

using namespace std;


//  Stores a histogram of kmer count values.


class merylHistogram {
public:
  merylHistogram();
  ~merylHistogram();

  void      addValue(kmvalu value) {

    if (value == 0)
      return;

    if (value == 1)
      _numUnique++;

    _numDistinct += 1;
    _numTotal    += value;

    if (value < _histMax)
      _hist[value]++;
    else
      _histBig[value]++;
  };

  void      clear(void);

  void      dump(stuffedBits *bits);
  void      dump(FILE        *outFile);

  void      load_v01(stuffedBits *bits);
  void      load_v03(stuffedBits *bits);
  void      load(stuffedBits *bits,    uint32 version);
  void      load(FILE        *inFile,  uint32 version);

  uint64    numUnique(void)                       { return(_numUnique);   };
  uint64    numDistinct(void)                     { return(_numDistinct); };
  uint64    numTotal(void)                        { return(_numTotal);    };

  uint32    histogramLength(void)                 { return(_histLen);     };
  uint64    histogramValue(uint32 i)              { return(_histVs[i]);   };
  uint64    histogramOccurrences(uint32 i)        { return(_histOs[i]);   };

private:
  uint64              _numUnique;
  uint64              _numDistinct;
  uint64              _numTotal;

  uint32              _histMax;    //  Max value that can be stored in _hist.
  uint64             *_hist;
  map<uint64, uint64> _histBig;    //  Values bigger than _histMax; <value,occurrances>

  uint64              _histLen;    //  If loaded from disk, this is the unpacked histogram.
  uint64             *_histVs;     //  The value this histogram entry is counting.
  uint64             *_histOs;     //  The number of occurrences of that value.
};


#endif  //  LIBKMER_HISTOGRAM_H
