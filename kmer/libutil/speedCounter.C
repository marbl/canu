#include <stdio.h>
#include <stdlib.h>

#include "util++.H"

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
