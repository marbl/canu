#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"
#include "bri++.H"




//  Command line options.  Only tmpPath needs to be global, and it can
//  be easily localized.
//
bool    beVerbose   = false;
u64bit  memoryLimit = 128 * 1024 * 1024;
char   *tmpPath     = 0L;



class aHitReader {
public:

  //  Open the file for reading, testing if it's binary or ascii input
  //
  aHitReader(char *filename) {
    errno = 0;
    theFile = fopen(filename, "r");
    if (theFile == 0L) {
      fprintf(stderr, "sortHits-- ERROR opening '%s': %s\n", filename, strerror(errno));
      exit(1);
    }

    char x = (char)fgetc(theFile);
    ungetc(x, theFile);

    isBinary = (x != '-');

    if (!isBinary)
      buffer = new char [1024];
  };

  ~aHitReader() {
    fclose(theFile);
    delete [] buffer;
  };

  bool   readHit(aHit &hit) {
    if (isBinary) {
      ahit_readBinary(&hit, theFile);
    } else {
      fgets(buffer, 1024, theFile);
      ahit_parseString(&hit, buffer);
    }

    return(feof(theFile) == false);
  };
private:
  FILE *theFile;
  char *buffer;
  bool  isBinary;
};



//  Write a bunch of hits to a temporary file (unlink the file after
//  it's opened) then allow those hits to be read back in.  Doesn't
//  need the aHitReader, as we use just the binary format.
//
class aHitTemporary {
public:
  aHitTemporary(aHit *hits, u32bit hitsLen) {
    theFile = makeTempFile(tmpPath);

    //  XXX:  Known bug on Tru64: fwrite() of data blocks > 2GB is broken

    u32bit  outputPos = 0;
    u32bit  outputLen = 1024 * 1024 / sizeof(aHit);

    while (outputPos < hitsLen) {
      errno = 0;
      outputPos += fwrite(hits, sizeof(aHit), hitsLen, theFile);
      if (errno) {
        fprintf(stderr, "ERROR: sortHits()-- Failed to write temporary file: %s\n", strerror(errno));
        exit(1);
      }

      //  XXX:  do we write one too many?

      if (outputPos + outputLen > hitsLen)
        outputLen = hitsLen - outputPos;
    }

    rewind(theFile);

    hit._direction = u32bitZERO;
    hit._qsIdx     = u32bitZERO;
    hit._dsIdx     = u32bitZERO;
    hit._dsLo      = u32bitZERO;
    hit._dsHi      = u32bitZERO;
    hit._covered   = u32bitZERO;
    hit._matched   = u32bitZERO;
    hit._numMers   = u32bitZERO;

    nextHit();
  };

  ~aHitTemporary() {
    fclose(theFile);
  };

  aHit  *theHit(void) {
    return(&hit);
  };
  void   nextHit(void) {
    if (hit._direction != ~u32bitZERO) {
      errno = 0;
      fread(&hit, sizeof(aHit), 1, theFile);
      if (errno) {
        fprintf(stderr, "ERROR: sortHits()-- Failed to read a hit: %s\n", strerror(errno));
        exit(1);
      }

      //  If we hit eof, this hit is invalid, and so are all future ones.  Set
      //  hit to be junk.
      //
      if (feof(theFile)) {
        hit._direction = ~u32bitZERO;
        hit._qsIdx     = ~u32bitZERO;
        hit._dsIdx     = ~u32bitZERO;
        hit._dsLo      = ~u32bitZERO;
        hit._dsHi      = ~u32bitZERO;
        hit._covered   = ~u32bitZERO;
        hit._matched   = ~u32bitZERO;
        hit._numMers   = ~u32bitZERO;
      }
    }
  };
private:
  FILE *theFile;
  aHit  hit;
};






int
hitcmp(const void *a, const void *b) {
  aHit  *A = (aHit *)a;
  aHit  *B = (aHit *)b;

  if (A->_dsIdx < B->_dsIdx) return(-1);
  if (A->_dsIdx > B->_dsIdx) return(1);
  if (A->_qsIdx < B->_qsIdx) return(-1);
  if (A->_qsIdx > B->_qsIdx) return(1);
  if (A->_dsLo  < B->_dsLo)  return(-1);
  if (A->_dsLo  > B->_dsLo)  return(1);
  return(0);
}




int
main(int argc, char **argv) {

  if (argc < 4) {
    fprintf(stderr, "usage: %s [-v] [-m memorylimit] [-t temppath] hitfile1 hitfile2 ... > sorted-hits\n", argv[0]);
    fprintf(stderr, "       memory limit is MB\n");
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-m", 2) == 0) {
      arg++;
      memoryLimit   = atoi(argv[arg]);
      memoryLimit <<= 20;
    } else if (strncmp(argv[arg], "-t", 2) == 0) {
      arg++;
      tmpPath = argv[arg];
    } else {
      //  Must be at the first file name.  Break.
      break;
    }
    arg++;
  }

  //  Allocate a bunch of spaces to store hits.
  //
  u64bit   hitsMax = memoryLimit / sizeof(aHit);
  u32bit   hitsPos = 0;
  aHit  *hits    = new aHit [hitsMax];

  u32bit            tmpFlen = 0;
  u32bit            tmpFmax = 1024;
  aHitTemporary   **tmpF    = new aHitTemporary * [tmpFmax];

  while (arg < argc) {
    aHitReader  *R = new aHitReader(argv[arg]);
    arg++;

    //  Read hits until we exhaust out space, then sort and dump to disk.
    //
    while (R->readHit(hits[hitsPos])) {
      hitsPos++;

      if (hitsPos == hitsMax) {
        qsort(hits, hitsPos, sizeof(aHit), hitcmp);

        if (tmpFlen >= tmpFmax) {
          tmpFmax *= 2;
          aHitTemporary **tmp = new aHitTemporary * [tmpFmax];
          memcpy(tmp, tmpF, sizeof(aHitTemporary) * tmpFlen);
          delete [] tmpF;
          tmpF = tmp;
        }

        tmpF[tmpFlen] = new aHitTemporary(hits, hitsPos);

        tmpFlen++;
        hitsPos = 0;
      }
    }

    delete R;
  }

  //  All done reading.  If we have stuff to sort, sort it.
  //
  if (hitsPos > 0)
    qsort(hits, hitsPos, sizeof(aHit), hitcmp);

  //  No temporary files?  Just write the hits and exit.  We're done.
  //
  if (tmpFlen == 0) {
    for (u32bit i=0; i<hitsPos; i++)
      ahit_printASCII(hits+i, stdout);
    exit(0);
  }


  //  We have at least one temporary file already on disk, so to make things
  //  easier, we write out the current set of hits and do an all disk merge.


  if (tmpFlen >= tmpFmax) {
    tmpFmax *= 2;
    aHitTemporary **tmp = new aHitTemporary * [tmpFmax];
    memcpy(tmp, tmpF, sizeof(aHitTemporary) * tmpFlen);
    delete [] tmpF;
    tmpF = tmp;
  }

  tmpF[tmpFlen] = new aHitTemporary(hits, hitsPos);

  tmpFlen++;



  //  While there is still input, merge to the output
  //
  bool moreInput = true;

  while (moreInput) {

    //  Pick the smallest hit -- if file [i] is finished, then hit[i]
    //  is bogus and all the values are set to maximal values.
    //
    u32bit smallestHit = 0;

    for (u32bit nh = smallestHit+1; nh < tmpFlen; nh++) {
      if (hitcmp(tmpF[smallestHit]->theHit(), tmpF[nh]->theHit()) > 0)
        smallestHit = nh;
    }

    //  If the smallest hit is invalid, we're done.  Otherwise, write
    //  the hit, and read a new one.
    //
    //
    if (tmpF[smallestHit]->theHit()->_direction == ~u32bitZERO) {
      moreInput = false;
    } else {
      ahit_printASCII(tmpF[smallestHit]->theHit(), stdout);
      tmpF[smallestHit]->nextHit();
    }
  }

  //  Should clean up, I know.

  return(0);
}
