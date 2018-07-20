
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
 *    kmer/libutil/speedCounter.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2006-OCT-22 to 2014-APR-11
 *      are Copyright 2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "speedCounter.H"

const char*
speedCounter::_spinr[4] = { "[|]", "[/]", "[-]", "[\\]" };

const char*
speedCounter::_liner[19] = { "[-         ]",
                             "[--        ]",
                             "[ --       ]",
                             "[  --      ]",
                             "[   --     ]",
                             "[    --    ]",
                             "[     --   ]",
                             "[      --  ]",
                             "[       -- ]",
                             "[        --]",
                             "[         -]",
                             "[        --]",
                             "[       -- ]",
                             "[      --  ]",
                             "[     --   ]",
                             "[    --    ]",
                             "[   --     ]",
                             "[  --      ]",
                             "[ --       ]" };


speedCounter::speedCounter(char const   *fmt,
                           double        unit,
                           uint64        freq,
                           bool          enabled) {
  _count     = 0;
  _draws     = 0;
  _unit      = unit;
  _freq      = freq;
  _startTime = getTime();
  _fmt       = fmt;
  _spin      = false;
  _line      = false;
  _enabled   = enabled;

  //  We use _draws instead of shifting _count just because it's
  //  simpler, and both methods need another variable anyway.

  //  Set all the bits below the hightest set in _freq --
  //  this allows us to do a super-fast test in tick().
  //
  _freq |= _freq >> 1;
  _freq |= _freq >> 2;
  _freq |= _freq >> 4;
  _freq |= _freq >> 8;
  _freq |= _freq >> 16;
  _freq |= _freq >> 32;
}

speedCounter::~speedCounter() {
  finish();
}
