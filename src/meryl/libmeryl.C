
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
 *    kmer/libmeryl/libmeryl.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-SEP-08 to 2004-APR-08
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-APR-30
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-23 to 2014-APR-11
 *      are Copyright 2005-2008,2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-JUL-01
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "libmeryl.H"

#include "AS_UTL_fileIO.H"
#include "AS_UTL_alloc.H"


//  Version 3 ??
//  Version 4 removed _histogramHuge, dynamically sizing it on write.

//                      0123456789012345
static char *ImagicV = "merylStreamIv04\n";
static char *ImagicX = "merylStreamIvXX\n";
static char *DmagicV = "merylStreamDv04\n";
static char *DmagicX = "merylStreamDvXX\n";
static char *PmagicV = "merylStreamPv04\n";
static char *PmagicX = "merylStreamPvXX\n";

merylStreamReader::merylStreamReader(const char *fn_, uint32 ms_) {

  if (fn_ == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  memset(_filename, 0, sizeof(char) * FILENAME_MAX);
  strcpy(_filename, fn_);

  //  Open the files
  //
  char *inpath = new char [strlen(_filename) + 8];

  sprintf(inpath, "%s.mcidx", _filename);
  _IDX = new bitPackedFile(inpath);

  sprintf(inpath, "%s.mcdat", _filename);
  _DAT = new bitPackedFile(inpath);

  sprintf(inpath, "%s.mcpos", _filename);
  if (AS_UTL_fileExists(inpath))
    _POS = new bitPackedFile(inpath);
  else
    _POS = 0L;

  delete [] inpath;

  //  Verify that they are what they should be, and read in the header
  //
  char    Imagic[16] = {0};
  char    Dmagic[16] = {0};
  char    Pmagic[16] = {0};
  bool    fail       = false;

  for (uint32 i=0; i<16; i++) {
    Imagic[i] = _IDX->getBits(8);
    Dmagic[i] = _DAT->getBits(8);
    if (_POS)
      Pmagic[i] = _POS->getBits(8);
  }

  if (strncmp(Imagic, ImagicX, 16) == 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx is an INCOMPLETE merylStream index file!\n", _filename);
    fail = true;
  }
  if (strncmp(Imagic, ImagicX, 13) != 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx is not a merylStream index file!\n", _filename);
    fail = true;
  }

  if (strncmp(Dmagic, DmagicX, 16) == 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcdat is an INCOMPLETE merylStream data file!\n", _filename);
    fail = true;
  }
  if (strncmp(Dmagic, DmagicX, 13) != 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcdat is not a merylStream data file!\n", _filename);
    fail = true;
  }

  if ((Imagic[13] != Dmagic[13]) ||
      (Imagic[14] != Dmagic[14])) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx and %s.mcdat are different versions!\n", _filename, _filename);
    fail = true;
  }

  if (_POS) {
    if (strncmp(Pmagic, PmagicX, 16) == 0) {
      fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcpos is an INCOMPLETE merylStream data file!\n", _filename);
      fail = true;
    }
    if (strncmp(Pmagic, PmagicX, 13) != 0) {
      fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcpos is not a merylStream data file!\n", _filename);
      fail = true;
    }
  }

  if (fail)
    exit(1);

  _idxIsPacked    = _IDX->getBits(32);
  _datIsPacked    = _IDX->getBits(32);
  _posIsPacked    = _IDX->getBits(32);

  _merSizeInBits  = _IDX->getBits(32) << 1;
  _merCompression = _IDX->getBits(32);
  _prefixSize     = _IDX->getBits(32);
  _merDataSize    = _merSizeInBits - _prefixSize;

  _numUnique      = _IDX->getBits(64);
  _numDistinct    = _IDX->getBits(64);
  _numTotal       = _IDX->getBits(64);

  _histogramPos      = 0;
  _histogramLen      = 0;
  _histogramMaxValue = 0;
  _histogram         = 0L;

  uint32 version = atoi(Imagic + 13);

  //  Versions earlier than four used a fixed-size histogram, stored at the start
  //  of the index.

  if (version < 4) {
    _histogramPos      = _IDX->tell();
    _histogramLen      = _IDX->getBits(64);  //  Previous _histogramHuge, now unused
    _histogramLen      = _IDX->getBits(64);
    _histogramMaxValue = _IDX->getBits(64);
    _histogram         = new uint64 [_histogramLen];

    for (uint32 i=0; i<_histogramLen; i++)
      _histogram[i] = _IDX->getBits(64);
  }

  //  Version 4 switched to a dynamically sized histogram, stored at the end
  //  of the index.

  else {
    _histogramPos      = _IDX->getBits(64);
    _histogramLen      = _IDX->getBits(64);
    _histogramMaxValue = _IDX->getBits(64);
    _histogram         = new uint64 [_histogramLen];

    uint64  position = _IDX->tell();

    _IDX->seek(_histogramPos);

    for (uint32 i=0; i<_histogramLen; i++)
      _histogram[i] = _IDX->getBits(64);

    _IDX->seek(position);
  }


  _thisBucket     = uint64ZERO;
  _thisBucketSize = getIDXnumber();
  _numBuckets     = uint64ONE << _prefixSize;

  _thisMer.setMerSize(_merSizeInBits >> 1);
  _thisMer.clear();
  _thisMerCount   = uint64ZERO;

  _thisMerPositionsMax = 0;
  _thisMerPositions    = 0L;

  _validMer       = true;

#ifdef SHOW_VARIABLES
  fprintf(stderr, "_merSizeInBits  = "F_U32"\n", _merSizeInBits);
  fprintf(stderr, "_merCompression = "F_U32"\n", _merCompression);
  fprintf(stderr, "_prefixSize     = "F_U32"\n", _prefixSize);
  fprintf(stderr, "_merDataSize    = "F_U32"\n", _merDataSize);
  fprintf(stderr, "_numUnique      = "F_U64"\n", _numUnique);
  fprintf(stderr, "_numDistinct    = "F_U64"\n", _numDistinct);
  fprintf(stderr, "_numTotal       = "F_U64"\n", _numTotal);
  fprintf(stderr, "_thisBucket     = "F_U64"\n", _thisBucket);
  fprintf(stderr, "_thisBucketSize = "F_U64"\n", _thisBucketSize);
  fprintf(stderr, "_thisMerCount   = "F_U64"\n", _thisMerCount);
#endif

  if ((ms_ > 0) && (_merSizeInBits >> 1 != ms_)) {
    fprintf(stderr, "merylStreamReader()-- ERROR: User requested mersize "F_U32" but '%s' is mersize "F_U32"\n",
            ms_, _filename, _merSizeInBits >> 1);
    exit(1);
  }
}


merylStreamReader::~merylStreamReader() {
  delete _IDX;
  delete _DAT;
  delete _POS;
  delete [] _thisMerPositions;
  delete [] _histogram;
}



bool
merylStreamReader::nextMer(void) {

  //  Use a while here, so that we skip buckets that are empty
  //
  while ((_thisBucketSize == 0) && (_thisBucket < _numBuckets)) {
    _thisBucketSize = getIDXnumber();
    _thisBucket++;
  }

  if (_thisBucket >= _numBuckets)
    return(_validMer = false);

  //  Before you get rid of the clear() -- if, say, the list of mers
  //  is sorted and we can shift the mer to make space for the new
  //  stuff -- make sure that nobody is calling reverseComplement()!
  //
  _thisMer.clear();
  _thisMer.readFromBitPackedFile(_DAT, _merDataSize);
  _thisMer.setBits(_merDataSize, _prefixSize, _thisBucket);

  _thisMerCount = getDATnumber();

  _thisBucketSize--;

  if (_POS) {
    if (_thisMerPositionsMax < _thisMerCount) {
      delete [] _thisMerPositions;
      _thisMerPositionsMax = _thisMerCount + 1024;
      _thisMerPositions    = new uint32 [_thisMerPositionsMax];
    }
    for (uint32 i=0; i<_thisMerCount; i++) {
      _thisMerPositions[i] = _POS->getBits(32);
    }
  }

  return(true);
}






merylStreamWriter::merylStreamWriter(const char *fn_,
                                     uint32 merSize,
                                     uint32 merComp,
                                     uint32 prefixSize,
                                     bool   positionsEnabled) {

  memset(_filename, 0, sizeof(char) * FILENAME_MAX);
  strcpy(_filename, fn_);

  char *outpath = new char [FILENAME_MAX];

  sprintf(outpath, "%s.mcidx.creating", _filename);
  _IDX = new bitPackedFile(outpath, 0, true);

  sprintf(outpath, "%s.mcdat.creating", _filename);
  _DAT = new bitPackedFile(outpath, 0, true);

  if (positionsEnabled) {
    sprintf(outpath, "%s.mcpos.creating", _filename);
    _POS = new bitPackedFile(outpath, 0, true);
  } else {
    _POS = 0L;
  }

  delete [] outpath;

  _idxIsPacked    = 1;
  _datIsPacked    = 1;
  _posIsPacked    = 0;

  _merSizeInBits  = merSize * 2;
  _merCompression = merComp;
  _prefixSize     = prefixSize;
  _merDataSize    = _merSizeInBits - _prefixSize;

  _thisBucket     = uint64ZERO;
  _thisBucketSize = uint64ZERO;
  _numBuckets     = uint64ONE << _prefixSize;

  _numUnique      = uint64ZERO;
  _numDistinct    = uint64ZERO;
  _numTotal       = uint64ZERO;

  _histogramPos      = 0;
  _histogramLen      = 1024;
  _histogramMaxValue = 0;
  _histogram         = new uint64 [_histogramLen];

  for (uint32 i=0; i<_histogramLen; i++)
    _histogram[i] = 0;

  _thisMerIsBits  = false;
  _thisMerIskMer  = false;

  _thisMer.setMerSize(_merSizeInBits >> 1);
  _thisMer.clear();

  _thisMerPre     = uint64ZERO;
  _thisMerMer     = uint64ZERO;

  _thisMerPreSize = prefixSize;
  _thisMerMerSize = 2 * merSize - prefixSize;

  _thisMerCount   = uint64ZERO;

  //  Initialize the index file.

  for (uint32 i=0; i<16; i++)
    _IDX->putBits(ImagicX[i], 8);

  _IDX->putBits(_idxIsPacked, 32);
  _IDX->putBits(_datIsPacked, 32);
  _IDX->putBits(_posIsPacked, 32);

  _IDX->putBits(_merSizeInBits >> 1, 32);
  _IDX->putBits(_merCompression, 32);
  _IDX->putBits(_prefixSize,  32);
  _IDX->putBits(_numUnique,   64);
  _IDX->putBits(_numDistinct, 64);
  _IDX->putBits(_numTotal,    64);

  _IDX->putBits(0, 64);        //  Offset to the histogram
  _IDX->putBits(0, 64);        //  Length of the histogram data
  _IDX->putBits(0, 64);        //  Max value seen in the histogram

  //  Initialize the data file.

  for (uint32 i=0; i<16; i++)
    _DAT->putBits(DmagicX[i], 8);

  //  Initialize the positions file.

  if (_POS)
    for (uint32 i=0; i<16; i++)
      _POS->putBits(PmagicX[i], 8);
}


merylStreamWriter::~merylStreamWriter() {

  writeMer();

  //  Finish writing the buckets.

  while (_thisBucket < _numBuckets + 2) {
    setIDXnumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  //  Save the position of the histogram

  _histogramPos = _IDX->tell();

  //  And write the histogram

  for (uint32 i=0; i<=_histogramMaxValue; i++)
    _IDX->putBits(_histogram[i], 64);

  //  Seek back to the start and rewrite the magic numbers.

  _IDX->seek(0);

  for (uint32 i=0; i<16; i++)
    _IDX->putBits(ImagicV[i], 8);

  _IDX->putBits(_idxIsPacked, 32);
  _IDX->putBits(_datIsPacked, 32);
  _IDX->putBits(_posIsPacked, 32);

  _IDX->putBits(_merSizeInBits >> 1, 32);
  _IDX->putBits(_merCompression, 32);
  _IDX->putBits(_prefixSize,  32);
  _IDX->putBits(_numUnique,   64);
  _IDX->putBits(_numDistinct, 64);
  _IDX->putBits(_numTotal,    64);

  _IDX->putBits(_histogramPos,        64);
  _IDX->putBits(_histogramMaxValue+1, 64);  //  The length of the data (includes 0)
  _IDX->putBits(_histogramMaxValue,   64);  //  The maximum value of the data

  delete    _IDX;
  delete [] _histogram;

  //  Seek back to the start of the data and rewrite the magic numbers.

  _DAT->seek(0);

  for (uint32 i=0; i<16; i++)
    _DAT->putBits(DmagicV[i], 8);

  delete _DAT;

  //  Seek back to the start of the positions and rewrite the magic numbers.

  if (_POS) {
    _POS->seek(0);

    for (uint32 i=0; i<16; i++)
      _POS->putBits(PmagicV[i], 8);
  }

  delete _POS;

  //  All done!  Rename our temporary outputs to final outputs.

  char *outpath = new char [FILENAME_MAX];
  char *finpath = new char [FILENAME_MAX];

  sprintf(outpath, "%s.mcidx.creating", _filename);
  sprintf(finpath, "%s.mcidx", _filename);
  rename(outpath, finpath);

  sprintf(outpath, "%s.mcdat.creating", _filename);
  sprintf(finpath, "%s.mcdat", _filename);
  rename(outpath, finpath);

  if (_POS) {
    sprintf(outpath, "%s.mcpos.creating", _filename);
    sprintf(finpath, "%s.mcpos", _filename);
    rename(outpath, finpath);
  }

  delete [] finpath;
  delete [] outpath;
}


void
merylStreamWriter::writeMer(void) {

  if (_thisMerCount == 0)
    return;

  _numTotal += _thisMerCount;
  _numDistinct++;

  if (_thisMerCount >= _histogramLen)
    resizeArray(_histogram, _histogramMaxValue, _histogramLen, _thisMerCount + 16384);

  _histogram[_thisMerCount]++;

  if (_histogramMaxValue < _thisMerCount)
    _histogramMaxValue = _thisMerCount;

  assert((_thisMerIsBits == false) || (_thisMerIskMer == false));

  if (_thisMerIsBits) {
    if (_thisMerCount == 1) {
      _DAT->putBits(_thisMerMer, _thisMerMerSize);
      setDATnumber(1);
      _thisBucketSize++;
      _numUnique++;
    } else {
      _DAT->putBits(_thisMerMer, _thisMerMerSize);
      setDATnumber(_thisMerCount);
      _thisBucketSize++;
    }

  } else {
    if (_thisMerCount == 1) {
      _thisMer.writeToBitPackedFile(_DAT, _merDataSize);
      setDATnumber(1);
      _thisBucketSize++;
      _numUnique++;
    } else if (_thisMerCount > 1) {
      _thisMer.writeToBitPackedFile(_DAT, _merDataSize);
      setDATnumber(_thisMerCount);
      _thisBucketSize++;
    }
  }
}



void
merylStreamWriter::addMer(kMer &mer, uint32 count, uint32 *positions) {
  uint64  val;

  if (_thisMerIskMer == false) {
    _thisMerIskMer = true;
    assert(_thisMerIsBits == false);
  }

  //  Fail if we see a smaller mer than last time.
  //
  if (mer < _thisMer) {
    char str[1024];
    fprintf(stderr, "merylStreamWriter::addMer()-- ERROR: your mer stream isn't sorted increasingly!\n");
    fprintf(stderr, "merylStreamWriter::addMer()-- last: %s\n", _thisMer.merToString(str));
    fprintf(stderr, "merylStreamWriter::addMer()-- this: %s\n", mer.merToString(str));
    exit(1);
  }

  //  If there was a position given, write it.
  //
  if (positions && _POS)
    for (uint32 i=0; i<count; i++)
      _POS->putBits(positions[i], 32);

  //  If the new mer is the same as the last one just increase the
  //  count.
  //
  if (mer == _thisMer) {
    _thisMerCount += count;
    return;
  }

  //  Write thisMer to disk.  If the count is zero, we don't write
  //  anything.  The count is zero for the first mer (all A) unless we
  //  add that mer, and if the silly user gives us a mer with zero
  //  count.
  //
  writeMer();

  //  If the new mer is in a different bucket from the last mer, write
  //  out some bucket counts.  We need a while loop (opposed to just
  //  writing one bucket) because we aren't guaranteed that the mers
  //  are in adjacent buckets.
  //
  val = mer.startOfMer(_prefixSize);

  while (_thisBucket < val) {
    setIDXnumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  //  Remember the new mer for the next time
  //
  _thisMer      = mer;
  _thisMerCount = count;
}



void
merylStreamWriter::addMer(uint64  prefix,  uint32 prefixBits,
                          uint64  mer,     uint32 merBits,
                          uint32  count,
                          uint32 *UNUSED(positions)) {

  if (_thisMerIsBits == false) {
    _thisMerIsBits = true;
    assert(_thisMerIskMer == false);
  }

  assert(prefixBits           == _prefixSize);
  assert(prefixBits           == _thisMerPreSize);
  assert(merBits              == _thisMerMerSize);
  assert(prefixBits + merBits == _merSizeInBits);

  if (((prefix <  _thisMerPre)) ||
      ((prefix <= _thisMerPre) && (mer < _thisMerMer))) {
    assert(0);
  }

  if ((prefix == _thisMerPre) &&
      (mer    == _thisMerMer)) {
    _thisMerCount += count;
    return;
  }

  writeMer();

  while (_thisBucket < prefix) {
    setIDXnumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  _thisMerPre   = prefix;
  _thisMerMer   = mer;
  _thisMerCount = count;
}
