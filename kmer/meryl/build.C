#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libbri.H"
#include "mcDescription.H"
#include "outputMer.H"
#include "time.H"

static mcDescription   mcd;

u64bit     *_chck;
u64bit     *_hash;

//  Stuff for sorting the list of mers.
//
//  This used to be u32bit, as the check shouldn't be more than 32
//  bits but it could.
//
typedef u64bit heapbit;

////////////////////////////////////////////////////////////////////////////////
//
//  BUILD
//
//

void
createHashTable(char *inputFile,
                bool  doForward,
                bool  doReverse,
                bool  doCanonical,
                bool  beVerbose) {

  if (beVerbose)
#ifdef TRUE64BIT
    fprintf(stderr, " 0) Allocating %4luMB for counting table.\n",
            (mcd._tableSizeInEntries) >> 18);
#else
    fprintf(stderr, " 0) Allocating %4lluMB for counting table.\n",
            (mcd._tableSizeInEntries) >> 18);
#endif

  u32bit *_ctbl = new u32bit [ mcd._tableSizeInEntries ];
  for (u64bit i=mcd._tableSizeInEntries; i--; )
    _ctbl[i] = 0;

  if (beVerbose)
    fprintf(stderr, " 1) Counting mers in buckets.\n");

  //
  //  XXX:  Huge stalls accessing _ctbl.
  //
  //  I think we can get by with 24 bits of count here; the fragstore
  //  seems to have only 9 million mers in the first bucket.  Should
  //  probably check for overflow, though.
  //  
  merStream          M(mcd._merSizeInBases, inputFile);
  speedCounter       C("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  while (M.nextMer()) {
    if (doCanonical) {
      doForward = (M.theFMer() <= M.theRMer());
      doReverse = !doForward;
    }

    if (doForward) {
      _ctbl[ mcd.HASH(M.theFMer()) ]++;
      mcd._actualNumberOfMers++;
      C.tick();
    }

    if (doReverse) {
      _ctbl[ mcd.HASH(M.theRMer()) ]++;
      mcd._actualNumberOfMers++;
      C.tick();
    }

    //fprintf(stderr, "F=0x%016lx R=0x%016lx\n", M.theFMer(), M.theRMer());
  }

  if (beVerbose) {
#ifdef TRUE64BIT
    fprintf(stderr, "\n    Found %lu mers.\n", mcd._actualNumberOfMers);
#else
    fprintf(stderr, "\n    Found %llu mers.\n", mcd._actualNumberOfMers);
#endif
  }

  //
  //  Allocate a PACKED array for the hash table.  This needs to be
  //  packed only for mcd._actualNumberOfMers > 4 billion, really.
  //

  //  Determine how many bits we need to hold the value
  //  mcd._actualNumberOfMers.....then....
  //
  //  This is mcd._actualNumberOfMers+1 because we need to store the
  //  first position after the last mer.  That is, if there are two
  //  mers, we will store that the first mer is at position 0, the
  //  second mer is at position 1, and the end of the second mer is at
  //  position 2.
  //
  if (mcd._hashWidth > 0) {
    u32bit hashWidthCheck  = 1;
    while ((mcd._actualNumberOfMers+1) > (u64bitONE << hashWidthCheck))
      hashWidthCheck++;
    if (mcd._hashWidth < hashWidthCheck) {
      fprintf(stderr, "ERROR:  Forced hash width of %d is too small.  Must be at least %d for this input.\n",
              mcd._hashWidth, hashWidthCheck);
      exit(1);
    }
  } else {
    mcd._hashWidth  = 1;
    while ((mcd._actualNumberOfMers+1) > (u64bitONE << mcd._hashWidth))
      mcd._hashWidth++;
  }

  //  ....allocate a hash table that is that many bits wide.
  //
  if (beVerbose)
#ifdef TRUE64BIT
    fprintf(stderr, " 2) Allocating %4luMB for hash table (%u bits wide).\n",
            ((mcd._tableSizeInEntries+1) * mcd._hashWidth / 64 + 2) >> 17, mcd._hashWidth);
#else
    fprintf(stderr, " 2) Allocating %4lluMB for hash table (%lu bits wide).\n",
            ((mcd._tableSizeInEntries+1) * mcd._hashWidth / 64 + 2) >> 17, mcd._hashWidth);
#endif
  _hash = new u64bit [(mcd._tableSizeInEntries+1) * mcd._hashWidth / 64 + 2];

  //
  //  Create the hash index using the counts.  The hash points
  //  to the end of the bucket; when we add a word, we move the
  //  hash bucket pointer down one.
  //
  //  When done, we can deallocate the counting table.
  //
  if (beVerbose)
    fprintf(stderr, " 3) Initializing the hash table.\n");

  u64bit i=0;
  u64bit j=0;
  u64bit c=0;

  while (i < mcd._tableSizeInEntries) {
    c += _ctbl[i++];
    setDecodedValue(_hash, j, mcd._hashWidth, c);
    j += mcd._hashWidth;
  }

  //  Add the location of the end of the table.  This is not
  //  modified when adding words, but is used to determine
  //  the size of the last bucket.
  //
  setDecodedValue(_hash, j, mcd._hashWidth, c);

  if (beVerbose)
    fprintf(stderr, " 4) Releasing counting table.\n");
  delete [] _ctbl;
}


void
verifyHashTable(void) {
  u64bit i=0, j=0, c=0, d=0;

  //fprintf(stderr, "    Verifying hash table.\n");

  while (i <= mcd._tableSizeInEntries) {
    c = getDecodedValue(_hash, j, mcd._hashWidth);

    if (c < d)
#ifdef TRUE64BIT
      fprintf(stderr, "ERROR:  Table[%lu] out of order.\n", i);
#else
      fprintf(stderr, "ERROR:  Table[%llu] out of order.\n", i);
#endif

    d = c;

    j += mcd._hashWidth;
    i++;
  }

  //fprintf(stderr, "    Verify finished.\n");
}


void
fillCheckTable(char *inputFile,
               bool  doForward,
               bool  doReverse,
               bool  doCanonical,
               bool  beVerbose) {

  //  Allocate space for mcd._actualNumberOfMers mers in the _chck array.
  //  This doesn't need to be cleared.
  //
  if (beVerbose)
#ifdef TRUE64BIT
    fprintf(stderr, " 4) Allocating %4luMB for check table (%u bits wide).\n",
            (mcd._actualNumberOfMers * mcd._chckBits / 64 + 1) >> 17, mcd._chckBits);
#else
    fprintf(stderr, " 4) Allocating %4lluMB for check table (%lu bits wide).\n",
            (mcd._actualNumberOfMers * mcd._chckBits / 64 + 1) >> 17, mcd._chckBits);
#endif
  _chck = new u64bit [mcd._actualNumberOfMers * mcd._chckBits / 64 + 1];

  if (beVerbose)
    fprintf(stderr, " 5) Filling mers into list.\n");

  merStream          M(mcd._merSizeInBases, inputFile);
  speedCounter       C("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  bool               moreMers = M.nextMer();

  while (moreMers) {
    u64bit fmer      = M.theFMer();
    u64bit rmer      = M.theRMer();

    u64bit  b;
    u64bit  cfp = 0;
    u64bit  crp = 0;

    if (doCanonical) {
      doForward = (fmer <= rmer);
      doReverse = !doForward;
    }

    if (doForward) {
      b = mcd.HASH(fmer) * mcd._hashWidth;
      cfp = preDecrementDecodedValue(_hash, b, mcd._hashWidth);
      C.tick();
    }

    if (doReverse) {
      b = mcd.HASH(rmer) * mcd._hashWidth;
      crp = preDecrementDecodedValue(_hash, b, mcd._hashWidth);
      C.tick();
    }

    //  this is here so that we have something to do while we are stalling on
    //  the cfp/crp stuff above
    //
    moreMers = M.nextMer();

    if (doForward) {
      cfp *= mcd._chckBits;
      setDecodedValue(_chck, cfp, mcd._chckBits, fmer & mcd._chckMask);
    }

    if (doReverse) {
      crp *= mcd._chckBits;
      setDecodedValue(_chck, crp, mcd._chckBits, rmer & mcd._chckMask);
    }
  }

  if (beVerbose)
    fprintf(stderr, "\n");
}






////////////////////////////////////////////////////////////////////////////////
//
//  OUTPUT
//
void
adjustHeap(heapbit *M, s64bit i, s64bit n) {
  heapbit   m = M[i];
  s64bit    j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && M[j] < M[j+1])
      j++;                   //  j is the larger child

    if (m >= M[j])           //  a position for M[i] has been found
      break;

    M[(j-1)/2] = M[j];       //  Move larger child up a level

    j = (j << 1) + 1;
  }

  M[(j-1)/2] = m;
}



void
sortAndOutput(char   *outfilename,
              u32bit  lowCount,
              u32bit  highCount,
              bool    beVerbose) {
  u64bit m     = u64bitONE << mcd._tableSizeInBits;
  u32bit count = 0;
  u32bit items = 0;

  heapbit *_sortedList    = 0L;
  u32bit   _sortedListMax = 0;
  u32bit   _sortedListLen = 0;


  if (beVerbose)
    fprintf(stderr, " 6) Writing output.\n");

  //fprintf(stderr, "lowCount  = %u\n", lowCount);
  //fprintf(stderr, "highCount = %u\n", highCount);

  //  Open the output files
  //
  char *outpath = new char [strlen(outfilename) + 17];

  sprintf(outpath, "%s.mcidx", outfilename);
  bitPackedFileWriter *IDX = new bitPackedFileWriter(outpath);

  sprintf(outpath, "%s.mcdat", outfilename);
  bitPackedFileWriter *DAT = new bitPackedFileWriter(outpath);

  delete [] outpath;


  //  Write the parameters to the DAT file.  It probably doesn't
  //  matter which file we write these to, as we can't really do
  //  random access anyway.
  //
  mcd.write(DAT);

  speedCounter  C(" %7.2f Mbuckets -- %5.2f Mbuckets/second\r", 1000000.0, 0x1fffff, beVerbose);

  //  For each bucket, sort it.  The output is done
  //  in the sort.
  //
  for (u64bit B=0, b=0; b<m; b++) {
    C.tick();

    u64bit st = getDecodedValue(_hash, B, mcd._hashWidth);
    B        += mcd._hashWidth;
    u64bit ed = getDecodedValue(_hash, B, mcd._hashWidth);

    if (ed < st)
#ifdef TRUE64BIT
      fprintf(stderr, "ERROR: Bucket %10lu starts at %10lu ends at %10lu\n", b, st, ed);
#else
      fprintf(stderr, "ERROR: Bucket %10llu starts at %10llu ends at %10llu\n", b, st, ed);
#endif

    _sortedListLen = (u32bit)(ed - st);

    count = 0;
    items = 0;

    if (_sortedListLen > 0) {

      //  Allocate more space, if we need to.
      //
      if (_sortedListLen > _sortedListMax) {
        delete [] _sortedList;
        _sortedList    = new heapbit [_sortedListLen + 1];
        _sortedListMax = _sortedListLen;
      }

      //  Unpack the check values
      //
      for (u64bit i=st, J=st*mcd._chckBits; i<ed; i++, J += mcd._chckBits)
        _sortedList[i-st] = getDecodedValue(_chck, J, mcd._chckBits);

      //  Sort if there is more than one item
      //
      if (_sortedListLen > 1) {

        //  Create the heap of lines.
        //
        for (s64bit t=(_sortedListLen-2)/2; t>=0; t--)
          adjustHeap(_sortedList, t, _sortedListLen);

        //  Interchange the new maximum with the element at the end of the tree
        //
        for (s64bit t=_sortedListLen-1; t>0; t--) {
          heapbit          tv = _sortedList[t];
          _sortedList[t]      = _sortedList[0];
          _sortedList[0]      = tv;

          adjustHeap(_sortedList, 0, t);
        }
      }


      //  Scan the list of sorted mers, counting them.  Whenever we 
      //  know the count, output it.
      //
      count = 1;
      if (_sortedListLen > 0) {
        for (u32bit t=1; t<_sortedListLen; t++) {
          if (_sortedList[t] != _sortedList[t-1]) {
            if ((lowCount <= count) && (count <= highCount)) {
              outputMer(DAT, mcd, b, _sortedList[t-1], count);
              items++;
            }
            count = 0;
          }

          count++;
        }

        if ((lowCount <= count) && (count <= highCount)) {
          outputMer(DAT, mcd, b, _sortedList[_sortedListLen-1], count);
          items++;
        }
      }
    }

    //  Output the index
    //
    IDX->putBits(items, 32);
  }

  delete DAT;
  delete IDX;

  if (beVerbose)
    fprintf(stderr, "\n");
}






void
build(char   *inputFile,
      char   *outputFile,
      u32bit  merSize,
      u32bit  tblSize,
      u32bit  hashSize,
      u32bit  lowCount,
      u32bit  highCount,
      bool    doForward,
      bool    doReverse,
      bool    doCanonical,
      bool    beVerbose) {

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no input file specified.\n");
    exit(1);
  }

  if (outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
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

  if (lowCount > highCount) {
    fprintf(stderr, "ERROR - lowCount > highCount??\n");
    exit(1);
  }

  mcd._merSizeInBases      = merSize;
  mcd._merSizeInBits       = mcd._merSizeInBases << 1;

  //  30jan03
  //  It would appear that we need at least 2 bits in the table.
  //
  if (mcd._merSizeInBits < tblSize + 2) {
    fprintf(stderr, "WARNING:  table is too big for the mer (mer is %u bits, table is %u bits).\n",
            mcd._merSizeInBits, tblSize);
    fprintf(stderr, "WARNING:  adjusting to table size (-t) of %u\n",
            mcd._merSizeInBits - 2);

    tblSize = mcd._merSizeInBits - 2;
  }

  mcd._tableSizeInBits     = tblSize;
  mcd._tableSizeInEntries  = u32bitONE << mcd._tableSizeInBits;

  mcd._chckBits            = mcd._merSizeInBits - mcd._tableSizeInBits;
  _chck                    = 0L;
  mcd._chckMask            = u64bitMASK(mcd._chckBits);

  mcd._hashWidth           = hashSize;
  _hash                    = 0L;
  mcd._hashMask            = u64bitMASK(mcd._tableSizeInBits);  //  unused?

  mcd._actualNumberOfMers  = 0;

  if (beVerbose)
    mcd.print(stderr);

  createHashTable(inputFile, doForward, doReverse, doCanonical, beVerbose);
  verifyHashTable();
  fillCheckTable(inputFile, doForward, doReverse, doCanonical, beVerbose);
  sortAndOutput(outputFile, lowCount, highCount, beVerbose);

  delete [] _chck;
  delete [] _hash;
}

