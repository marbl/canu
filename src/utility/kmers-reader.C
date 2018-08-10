
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#include "kmers.H"
#include "bits.H"

#include "files.H"




kmerCountFileReader::kmerCountFileReader(const char *inputName,
                                         bool        ignoreStats,
                                         bool        beVerbose) {
  char   N[FILENAME_MAX+1];

  //  Save the input name for later use, but fail if
  //  the index file isn't found.

  strncpy(_inName, inputName, FILENAME_MAX);

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _inName);

  if (fileExists(N) == false)
    fprintf(stderr, "ERROR: '%s' doesn't appear to be a meryl input; file '%s' doesn't exist.\n",
            _inName, N), exit(1);

  //  Open the master index and initialize from it.

  stuffedBits  *masterIndex = new stuffedBits(N);

  uint64  m1 = masterIndex->getBinary(64);
  uint64  m2 = masterIndex->getBinary(64);

  if ((m1 != 0x646e496c7972656dllu) ||  //  merylInd
      (m2 != 0x31302e765f5f7865llu))    //  ex__v.01
    fprintf(stderr, "ERROR: '%s' doesn't look like a meryl input; file '%s' fails magic number check.\n",
            _inName, N), exit(1);

  _prefixSize    = masterIndex->getBinary(32);
  _suffixSize    = masterIndex->getBinary(32);

  _numFilesBits  = masterIndex->getBinary(32);
  _numBlocksBits = masterIndex->getBinary(32);

  _numFiles      = (uint64)1 << _numFilesBits;
  _numBlocks     = (uint64)1 << _numBlocksBits;

  _datFile       = NULL;

  _block         = new kmerCountFileReaderBlock();
  _blockIndex    = NULL;

  _kmer          = kmer();
  _count         = 0;

  _prefix        = 0;

  _activeMer     = 0;
  _activeFile    = 0;

  _nKmers        = 0;
  _nKmersMax     = 1024;
  _suffixes      = new uint64 [_nKmersMax];
  _counts        = new uint32 [_nKmersMax];

  if (ignoreStats == false)
    _stats.load(masterIndex);

  delete masterIndex;

  //  Check and setup the mer size if needed.

  uint32  merSize = (_prefixSize + _suffixSize) / 2;

  if (beVerbose) {
    char    m[17] = { 0 };

    for (uint32 i=0, s=0; i<8; i++, s+=8) {
      m[i + 0] = (m1 >> s) & 0xff;
      m[i + 8] = (m2 >> s) & 0xff;
    }

    fprintf(stderr, "Opened '%s'.\n", _inName);
    fprintf(stderr, "  magic          0x%016lx%016lx '%s'\n", m1, m2, m);
    fprintf(stderr, "  prefixSize     %u\n", _prefixSize);
    fprintf(stderr, "  suffixSize     %u\n", _suffixSize);
    fprintf(stderr, "  numFilesBits   %u (%u files)\n", _numFilesBits, _numFiles);
    fprintf(stderr, "  numBlocksBits  %u (%u blocks)\n", _numBlocksBits, _numBlocks);
  }

  if (kmer::merSize() == 0)    //  If the global kmer size isn't set yet, set it.
    kmer::setSize(merSize);    //  Then make sure all files are the same.

  if (kmer::merSize() != merSize)
    fprintf(stderr, "mer size mismatch, can't process this set of files.\n"), exit(1);
}



kmerCountFileReader::~kmerCountFileReader() {

  delete [] _blockIndex;

  delete [] _suffixes;
  delete [] _counts;

  AS_UTL_closeFile(_datFile);

  delete    _block;
}



void
kmerCountFileReader::loadBlockIndex(void) {

  if (_blockIndex != NULL)
    return;

  _blockIndex = new kmerCountFileIndex [_numFiles * _numBlocks];

  for (uint32 ii=0; ii<_numFiles; ii++) {
    char  *idxname = constructBlockName(_inName, ii, _numFiles, 0, true);
    FILE  *idxfile = AS_UTL_openInputFile(idxname);

    loadFromFile(_blockIndex + _numBlocks * ii, "kmerCountFileReader::blockIndex", _numBlocks, idxfile);

    AS_UTL_closeFile(idxfile, idxname);

    delete [] idxname;
  }
}



//  Like loadBlock, but just reports all blocks in the file, ignoring
//  the kmer data.
//
void
dumpMerylDataFile(char *name) {

  if (fileExists(name) == false)
    fprintf(stderr, "ERROR: '%s' doesn't exist.  Can't dump it.\n",
            name), exit(1);

  FILE                     *F = AS_UTL_openInputFile(name);
  stuffedBits              *D = new stuffedBits;

  fprintf(stdout, "            prefix   nKmers kCode uBits bBits                 k1 cCode                 c1                 c2\n");
  fprintf(stdout, "------------------ -------- ----- ----- ----- ------------------ ----- ------------------ ------------------\n");

  while (D->loadFromFile(F)) {
    uint64 position   = D->getPosition();

    uint64 m1         = D->getBinary(64);
    uint64 m2         = D->getBinary(64);

    uint64 prefix     = D->getBinary(64);
    uint64 nKmers     = D->getBinary(64);

    uint8  kCode      = D->getBinary(8);
    uint32 unaryBits  = D->getBinary(32);
    uint32 binaryBits = D->getBinary(32);
    uint64 k1         = D->getBinary(64);

    uint8  cCode      = D->getBinary(8);
    uint64 c1         = D->getBinary(64);
    uint64 c2         = D->getBinary(64);

    if ((m1 != 0x7461446c7972656dllu) ||
        (m2 != 0x0a3030656c694661llu)) {
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Magic number mismatch at position " F_U64 ".\n", position);
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x0a3030656c694661 got 0x%016" F_X64P "\n", m2);
      exit(1);
    }

    fprintf(stdout, "0x%016lx %8lu %5u %5u %5u 0x%016lx %5u 0x%016lx 0x%016lx\n",
            prefix, nKmers, kCode, unaryBits, binaryBits, k1, cCode, c1, c2);
  }

  delete D;

  AS_UTL_closeFile(F);
}



bool
kmerCountFileReader::nextMer(void) {

  _activeMer++;

  //  If we've still got data, just update and get outta here.
  //  Otherwise, we need to load another block.

  if (_activeMer < _nKmers) {
    _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
    _count = _counts[_activeMer];
    return(true);
  }

  //  Make sure all files are opened.

 loadAgain:
  if (_datFile == NULL)
    _datFile = openInputBlock(_inName, _activeFile, _numFiles);

  //  Load blocks.

  bool loaded = _block->loadBlock(_datFile, _activeFile);

  //  If nothing loaded. open a new file and try again.

  if (loaded == false) {
    AS_UTL_closeFile(_datFile);

    _activeFile++;

    if (_numFiles <= _activeFile)
      return(false);

    goto loadAgain;
  }

  //  Got a block!  Stash what we loaded.

  _prefix = _block->prefix();
  _nKmers = _block->nKmers();

#ifdef SHOW_LOAD
  fprintf(stdout, "LOADED prefix %016lx nKmers %lu\n", _prefix, _nKmers);
#endif

  //  Make sure we have space for the decoded data

  resizeArrayPair(_suffixes, _counts, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

  //  Decode the block into _OUR_ space.
  //
  //  decodeBlock() marks the block as having no data, so the next time we loadBlock() it will
  //  read more data from disk.  For blocks that don't get decoded, they retain whatever was
  //  loaded, and do not load another block in loadBlock().

  _block->decodeBlock(_suffixes, _counts);

  //  But if no kmers in this block, load another block.  Sadly, the block must always
  //  be decoded, otherwise, the load will not load a new block.

  if (_nKmers == 0)
    goto loadAgain;

  //  Reset iteration, and load the first kmer.

  _activeMer = 0;

  _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
  _count = _counts[_activeMer];

  return(true);
}
