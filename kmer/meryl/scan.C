#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bri++.H"
#include "meryl.H"



#error  Obsolete!  This file is going to be turned into merylScanner



"-S:  Given a table, and a list of mers, print the count for each mer.\n"
"        -s tblprefix  (use these tables)\n"
"        -q mers.fasta (get the list of mers from here; the query)"
"\n"
"        -o            (output file name)\n"
"\n"
"        -d            (include the defline in the output)\n"
"        -e            (state the mer explicitly)\n"
"\n"
"        -c            (basic output; for each mer, just the count)\n"
"        -p            (position output; for each mer, the position of the mer and the count)\n"
"        -a            (count for every base)\n"
"\n"
"     -d and -e only apply to -c and -p\n"
"     -a prints the count for every base, with deflines\n"
"\n"
"\n"




int
sortByMerHelper(const void *a, const void *b) {
  const mcMer *A = *((mcMer * const *)a);
  const mcMer *B = *((mcMer * const *)b);

  if (A->mer < B->mer)  return(-1);
  if (A->mer > B->mer)  return(1);
  return(0);
}

void
sortByMer(mcMer **m, u32bit l) {
  qsort(m, l, sizeof(mcMer *), sortByMerHelper);
}



int
sortByPositionHelper(const void *a, const void *b) {
  const mcMer *A = *((mcMer * const *)a);
  const mcMer *B = *((mcMer * const *)b);

  if (A->sequence < B->sequence)  return(-1);
  if (A->sequence > B->sequence)  return(1);
  if (A->position < B->position)  return(-1);
  if (A->position > B->position)  return(1);
  return(0);
}

void
sortByPosition(mcMer **m, u32bit l) {
  qsort(m, l, sizeof(mcMer*), sortByPositionHelper);
}



void
scan(merylArgs *args) {

  if (args->queryFile == 0L) {
    fprintf(stderr, "ERROR - no query file specified.\n");
    exit(1);
  }

  if (args->inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }

  if (!args->outputCount && !args->outputAll && !args->outputPosition) {
    fprintf(stderr, "ERROR:  You need to specify an output format: -c, -a or -p\n");
    exit(1);
  }



  merylStreamReader  *R = new merylStreamReader(args->inputFile);



  //  Open the output file
  //
  FILE *O = fopen(args->outputFile, "w");



  u32bit   _mersLen = 0;
  u32bit   _mersMax = 1048576;
  mcMer  **_mers    = new mcMer* [_mersMax];

  u32bit   _defLinesMax     = 0;
  char   **_defLines        = 0L;
  u32bit  *_seqLengths      = 0L;
  u32bit   _numberOfSeqs    = 0;

  if (args->includeDefLine || args->outputAll) {
    _defLinesMax     = 16384;
    _defLines        = new char*  [_defLinesMax];
    _seqLengths      = new u32bit [_defLinesMax];
    _numberOfSeqs    = 0;


    //  the zeroth sequence is not defined in FastAstream.
    //
    _defLines[0]   = 0L;
    _seqLengths[0] = 0;


    if (args->beVerbose)
      fprintf(stderr, " 0) Reading deflines and sequence lengths.\n");

    FastAstream    F(args->queryFile);
    unsigned char  ch = 255;

    while (ch != 0) {
      ch = F.nextSymbol();

      //  Got a new defline.  Save it.
      //
      if (ch == 254) {
        _numberOfSeqs++;

        //  Allocate more space?
        //
        if (F.theSequenceNumber() >= _defLinesMax) {
          _defLinesMax *= 2;

          char   **d = new char*  [_defLinesMax];
          u32bit  *l = new u32bit [_defLinesMax];

          for (u64bit z=F.theSequenceNumber(); z--; ) {
            d[z] = _defLines[z];
            l[z] = _seqLengths[z];
          }

          delete _defLines;
          delete _seqLengths;

          _defLines   = d;
          _seqLengths = l;
        }

        //fprintf(stderr, "%2lu] %s\n", F.theSequenceNumber(), F.theDefLine());

        //  Add the defline to our list
        //
        _defLines[F.theSequenceNumber()] = new char [strlen(F.theDefLine()) + 1];
        strcpy(_defLines[F.theSequenceNumber()], F.theDefLine());

        //  Clear the length
        //
        _seqLengths[F.theSequenceNumber()] = 0;
      } else {

        //  Increment the length (if it isn't end of file)
        //
        if ((ch > 0) && (ch < 254))
          _seqLengths[F.theSequenceNumber()]++;
      }
    }
  }

  if (args->beVerbose)
    fprintf(stderr, " 1) Reading mers.\n");

  //
  //  Read all the mers, sort them, then ask the counted sequence what
  //  the counts are.
  //
  merStream  *M = new merStream(mcd._merSizeInBases, args->queryFile);
  while (M->nextMer()) {

    //  Do we need more?
    //
    if (_mersLen >= _mersMax) {
      _mersMax <<= 1;

      mcMer **m = new mcMer* [_mersMax];

      for (u32bit i=_mersLen; i--; )
        m[i] = _mers[i];

      delete _mers;
      _mers = m;
    }

    _mers[_mersLen] = new mcMer;

    _mers[_mersLen]->position = (u32bit)M->thePosition();
    _mers[_mersLen]->sequence = (u32bit)M->theSequenceNumber();
    _mers[_mersLen]->count    = 0;

    if (args->doForward) {
      _mers[_mersLen]->mer = M->theFMer();
    } else if (args->doReverse) {
      _mers[_mersLen]->mer = M->theRMer();
    } else {
      _mers[_mersLen]->mer = M->theFMer();
      if (M->theRMer() < _mers[_mersLen]->mer)
        _mers[_mersLen]->mer = M->theRMer();
    }

    _mersLen++;
  }
  delete M;


  if (args->beVerbose) {
    fprintf(stderr, "    Found "u32bitFMT" mers.\n", _mersLen);
    fprintf(stderr, " 2) Sorting by mer.\n");
  }

  //  Sort mers.
  //
  sortByMer(_mers, _mersLen);


  if (args->beVerbose)
    fprintf(stderr, " 3) Scanning the mcFiles.\n");


  //
  //  Determine the count for each mer
  //
  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  for (u32bit i=0; i<_mersLen; i++) {
    if ((args->beVerbose) && ((i & 0xfff) == 0)) {
      fprintf(stderr, "%8.5f%% done, "u64bitFMT" bytes read.\r",
              100.0 * (double)i / (double)_mersLen,
              B->bitsRead() >> 3);
      fflush(stderr);
    }
    B->read(_mers[i]);
    B->scan(_mers[i]);
  }

  if (args->beVerbose)
    fprintf(stderr, "\n");

  if (args->beVerbose)
    fprintf(stderr, " 4) Sorting by position.\n");

  //  Sort again, putting things back into the original ordering
  //
  sortByPosition(_mers, _mersLen);


  //  Dump the output
  //
  if (args->beVerbose)
    fprintf(stderr, " 5) Writing the output.\n");

  u32bit  style = 0;
  u32bit  lastS = 0;

  if (args->includeDefLine)  style += 1;
  if (args->includeMer)      style += 2;

  if (args->outputCount) {
    switch (style) {
    case 0:
      for (u32bit i=0; i<_mersLen; i++) {
        fprintf(O, u32bitFMT"\n", _mers[i]->count);
      }
      break;
    case 1:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
        fprintf(O, u32bitFMT"\n", _mers[i]->count);
      }
      break;
    case 2:
      for (u32bit i=0; i<_mersLen; i++) {
        fprintf(O, u64bitHEX"\t"u32bitFMT"\n", _mers[i]->mer, _mers[i]->count);
      }
      break;
    case 3:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
        fprintf(O, u64bitHEX"\t"u32bitFMT"\n", _mers[i]->mer, _mers[i]->count);
      }
      break;
    }
  }



  if (args->outputAll) {
    u32bit thisMer = 0;  //  mer in the list
    u32bit thisPos = 0;  //  position in the current sequence
    u32bit thisSeq = 0;

    while (thisSeq<_numberOfSeqs) {

      //  Sequences are from 1 to n, inclusive, so we increment here.
      //
      thisSeq++;

      //  Always print the defline
      //
      fprintf(O, ">%s\n", _defLines[thisSeq]);

      //  If this sequence is missing, print all blank mers
      //
      if ((thisMer >= _mersLen) ||
          (thisSeq != _mers[thisMer]->sequence)) {
        for (thisPos=0; thisPos<_seqLengths[thisSeq]; thisPos++)
          fprintf(O, "-1\n");
      } else {
        //  Otherwise, print the sequence.
        //
        thisPos = 0;

        while ((thisPos  < _seqLengths[thisSeq]) &&
               (thisMer  < _mersLen) &&
               (thisSeq == _mers[thisMer]->sequence)) {

          //  Print any missing mers -- mers between
          //  thisPos and _mers[i].position
          //
          while (thisPos < _mers[thisMer]->position) {
            fprintf(O, "-1\n");
            thisPos++;
          }

          //  Print the mer
          //
          fprintf(O, u32bitFMT"\n", _mers[thisMer]->count);

          //  Next mer, next position
          //
          thisMer++;
          thisPos++;
        }

        //  We have finished printing all mers in the sequence.  Fill in any at the end.
        //
        if ((thisMer >= _mersLen) ||
            (thisSeq != _mers[thisMer]->sequence)) {
          for (; thisPos<_seqLengths[thisSeq]; thisPos++) {
            fprintf(O, "-1\n");
          }
        }
      }
    }
  }
  


  if (args->outputPosition) {
    switch (style) {
    case 0:
      for (u32bit i=0; i<_mersLen; i++) {
        fprintf(O, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
      }
      break;
    case 1:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
        fprintf(O, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
      }
      break;
    case 2:
      for (u32bit i=0; i<_mersLen; i++) {
        fprintf(O, u64bitHEX"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                _mers[i]->mer, _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
      }
      break;
    case 3:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
        fprintf(O, u64bitHEX"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                _mers[i]->mer, _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
      }
      break;
    }
  }




  //  Close stuff
  //
  fclose(O);

  delete DAT;
  delete IDX;
}
