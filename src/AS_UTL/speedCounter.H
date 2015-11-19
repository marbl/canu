
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
 *    kmer/libutil/speedCounter.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2003-MAY-06
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-12 to 2014-APR-11
 *      are Copyright 2005-2006,2012-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SPEEDCOUNTER_H
#define SPEEDCOUNTER_H

#include "AS_global.H"
#include "timeAndSize.H"

class speedCounter {
public:
  //  fmt specifies the status format.  An example:
  //    "    %8f [unit]things (%8.5f [unit]things/sec)\r"
  //
  speedCounter(char const   *fmt,
               double        unit,
               uint64        freq,
               bool          enabled=true);
  ~speedCounter();

  void   enableSpinner(void) { _spin = true; };
  void   enableLiner(void)   { _line = true; };

  bool   tick(void) {
    if (_enabled && ((++_count & _freq) == uint64ZERO)) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      _draws++;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
      return(true);
    }
    return(false);
  };

  bool   tick(uint64 increment) {
    if (_enabled == false)
      return(false);

    _count += increment;
    if ((_count & _freq) == uint64ZERO) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      _draws++;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
      return(true);
    }
    return(false);
  };

  void   finish(void) {
    if (_enabled && (_count >= _freq)) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fprintf(stderr, "\n");
      fflush(stderr);
    }
    _count = 0;
  };

private:
  static const char  *_spinr[4];
  static const char  *_liner[19];
  uint64              _count;
  uint64              _draws;
  double              _unit;
  uint64              _freq;
  double              _startTime;
  char const         *_fmt;
  bool                _spin;
  bool                _line;
  bool                _enabled;
};


#endif  //  SPEEDCOUNTER_H
