#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libbri.H"
#include "mcBucket.H"
#include "mcDescription.H"

static mcDescription mcd;

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
scan(char   *queryFile,
     char   *inputFile,
     char   *outputFile,
     bool    includeDefLine,
     bool    includeMer,
     bool    doForward,
     bool    doReverse,
     bool    doCanonical,
     bool    outputCount,
     bool    outputAll,
     bool    outputPosition,
     bool    beVerbose) {

  if (queryFile == 0L) {
    fprintf(stderr, "ERROR - no query file specified.\n");
    exit(1);
  }

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  if (outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }

  if (!outputCount && !outputAll && !outputPosition) {
    fprintf(stderr, "ERROR:  You need to specify an output format: -c, -a or -p\n");
    exit(1);
  }

  //  these should never happen, unles main() is broken.
  if ((doForward == false) && (doReverse == false) && (doCanonical == false)) {
    fprintf(stderr, "ERROR - need to specify at least one of -f, -r, -C\n");
    exit(1);
  }
  if ((doForward && doReverse) || (doForward && doCanonical) || (doReverse && doCanonical)) {
    fprintf(stderr, "ERROR - only one of -f, -r and -C may be specified!\n");
    exit(1);
  }



  //  Open the counted sequence files
  //
  char *inpath = new char [strlen(inputFile) + 17];

  sprintf(inpath, "%s.mcidx", inputFile);
  bitPackedFileReader *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", inputFile);
  bitPackedFileReader *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;


  //
  //  Read the parameters
  //
  mcd.read(DAT);


  //  Open the output file
  //
  FILE *O = fopen(outputFile, "w");



  u32bit   _mersLen = 0;
  u32bit   _mersMax = 1048576;
  mcMer  **_mers    = new mcMer* [_mersMax];

  u32bit   _defLinesMax     = 0;
  char   **_defLines        = 0L;
  u32bit  *_seqLengths      = 0L;
  u32bit   _numberOfSeqs    = 0;

  if (includeDefLine || outputAll) {
    _defLinesMax     = 16384;
    _defLines        = new char*  [_defLinesMax];
    _seqLengths      = new u32bit [_defLinesMax];
    _numberOfSeqs    = 0;


    //  the zeroth sequence is not defined in FastAstream.
    //
    _defLines[0]   = 0L;
    _seqLengths[0] = 0;


    if (beVerbose)
      fprintf(stderr, " 0) Reading deflines and sequence lengths.\n");

    FastAstream    F(queryFile);
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

  if (beVerbose)
    fprintf(stderr, " 1) Reading mers.\n");

  //
  //  Read all the mers, sort them, then ask the counted sequence what
  //  the counts are.
  //
  merStream  *M = new merStream(mcd._merSizeInBases, queryFile);
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

    if (doForward) {
      _mers[_mersLen]->mer = M->theFMer();
    } else if (doReverse) {
      _mers[_mersLen]->mer = M->theRMer();
    } else {
      _mers[_mersLen]->mer = M->theFMer();
      if (M->theRMer() < _mers[_mersLen]->mer)
        _mers[_mersLen]->mer = M->theRMer();
    }

    _mersLen++;
  }
  delete M;


  if (beVerbose) {
#ifdef TRUE64BIT
    fprintf(stderr, "    Found %u mers.\n", _mersLen);
#else
    fprintf(stderr, "    Found %lu mers.\n", _mersLen);
#endif
    fprintf(stderr, " 2) Sorting by mer.\n");
  }

  //  Sort mers.
  //
  sortByMer(_mers, _mersLen);


  if (beVerbose)
    fprintf(stderr, " 3) Scanning the mcFiles.\n");


  //
  //  Determine the count for each mer
  //
  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  for (u32bit i=0; i<_mersLen; i++) {
    if ((beVerbose) && ((i & 0xffff) == 0)) {
#ifdef TRUE64BIT
      fprintf(stderr, "%7u/%7u  %lu bits read\r", i, _mersLen, B->bitsRead());
#else
      fprintf(stderr, "%7lu/%7lu  %llu bits read\r", i, _mersLen, B->bitsRead());
#endif
      fflush(stderr);
    }
    B->read(_mers[i]);
    B->scan(_mers[i]);
  }

  if (beVerbose)
    fprintf(stderr, "\n");

  if (beVerbose)
    fprintf(stderr, " 4) Sorting by position.\n");

  //  Sort again, putting things back into the original ordering
  //
  sortByPosition(_mers, _mersLen);


  //  Dump the output
  //
  if (beVerbose)
    fprintf(stderr, " 5) Writing the output.\n");

  u32bit  style = 0;
  u32bit  lastS = 0;

  if (includeDefLine)  style += 1;
  if (includeMer)      style += 2;

  if (outputCount) {
    switch (style) {
    case 0:
      for (u32bit i=0; i<_mersLen; i++) {
#ifdef TRUE64BIT
        fprintf(O, "%u\n", _mers[i]->count);
#else
        fprintf(O, "%lu\n", _mers[i]->count);
#endif
      }
      break;
    case 1:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
#ifdef TRUE64BIT
        fprintf(O, "%u\n", _mers[i]->count);
#else
        fprintf(O, "%lu\n", _mers[i]->count);
#endif
      }
      break;
    case 2:
      for (u32bit i=0; i<_mersLen; i++) {
#ifdef TRUE64BIT
        fprintf(O, "0x%016lx\t%u\n", _mers[i]->mer, _mers[i]->count);
#else
        fprintf(O, "0x%016llx\t%lu\n", _mers[i]->mer, _mers[i]->count);
#endif
      }
      break;
    case 3:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
#ifdef TRUE64BIT
        fprintf(O, "0x%016lx\t%u\n", _mers[i]->mer, _mers[i]->count);
#else
        fprintf(O, "0x%016llx\t%lu\n", _mers[i]->mer, _mers[i]->count);
#endif
      }
      break;
    }
  }



  if (outputAll) {
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
#ifdef TRUE64BIT
          fprintf(O, "%u\n", _mers[thisMer]->count);
#else
          fprintf(O, "%lu\n", _mers[thisMer]->count);
#endif

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
  


  if (outputPosition) {
    switch (style) {
    case 0:
      for (u32bit i=0; i<_mersLen; i++) {
#ifdef TRUE64BIT
        fprintf(O, "%u\t%u\t%u\n",
                _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#else
        fprintf(O, "%lu\t%lu\t%lu\n",
                _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#endif
      }
      break;
    case 1:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
#ifdef TRUE64BIT
        fprintf(O, "%u\t%u\t%u\n",
                _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#else
        fprintf(O, "%lu\t%lu\t%lu\n",
                _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#endif
      }
      break;
    case 2:
      for (u32bit i=0; i<_mersLen; i++) {
#ifdef TRUE64BIT
        fprintf(O, "0x%016lx\t%u\t%u\t%u\n",
                _mers[i]->mer, _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#else
        fprintf(O, "0x%016llx\t%lu\t%lu\t%lu\n",
                _mers[i]->mer, _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#endif
      }
      break;
    case 3:
      for (u32bit i=0; i<_mersLen; i++) {
        while (lastS < _mers[i]->sequence) {
          lastS++;
          fprintf(O, ">%s\n", _defLines[lastS]);
        }
#ifdef TRUE64BIT
        fprintf(O, "0x%016lx\t%u\t%u\t%u\n",
                _mers[i]->mer, _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#else
        fprintf(O, "0x%016llx\t%lu\t%lu\t%lu\n",
                _mers[i]->mer, _mers[i]->sequence, _mers[i]->position, _mers[i]->count);
#endif
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
