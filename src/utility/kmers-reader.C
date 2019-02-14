
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

#include "kmers.H"
#include "bits.H"

#include "files.H"



//  Clear all members and allocate buffers.
void
kmerCountFileReader::initializeFromMasterI_v00(void) {

  _prefixSize    = 0;
  _suffixSize    = 0;

  _numFilesBits  = 0;
  _numBlocksBits = 0;

  _numFiles      = 0;
  _numBlocks     = 0;

  _stats         = NULL;

  _datFile       = NULL;

  _block         = new kmerCountFileReaderBlock();
  _blockIndex    = NULL;

  _kmer          = kmer();
  _value         = 0;

  _prefix        = 0;

  _activeMer     = 0;
  _activeFile    = 0;

  _threadFile    = UINT32_MAX;

  _nKmers        = 0;
  _nKmersMax     = 1024;
  _suffixes      = new uint64 [_nKmersMax];
  _values        = new uint64 [_nKmersMax];
}



//  Initialize for the original.
void
kmerCountFileReader::initializeFromMasterI_v01(stuffedBits  *masterIndex,
                                               bool          doInitialize) {

  if (doInitialize == true) {
    initializeFromMasterI_v00();

    _prefixSize    = masterIndex->getBinary(32);
    _suffixSize    = masterIndex->getBinary(32);

    _numFilesBits  = masterIndex->getBinary(32);
    _numBlocksBits = masterIndex->getBinary(32);

    _numFiles      = (uint64)1 << _numFilesBits;    //  The same for all formats, but
    _numBlocks     = (uint64)1 << _numBlocksBits;   //  awkward to do outside of here.
  }

  //  If we didn't initialize, set the file position to the start
  //  of the statistics.
  else {
    masterIndex->setPosition(64 + 64 + 32 + 32 + 32 + 32);
  }
}



//  Initialize for the format that includes multi sets.
void
kmerCountFileReader::initializeFromMasterI_v02(stuffedBits  *masterIndex,
                                               bool          doInitialize) {

  if (doInitialize == true) {
    initializeFromMasterI_v00();

    _prefixSize    = masterIndex->getBinary(32);
    _suffixSize    = masterIndex->getBinary(32);

    _numFilesBits  = masterIndex->getBinary(32);
    _numBlocksBits = masterIndex->getBinary(32);

    uint32 flags   = masterIndex->getBinary(32);

    _isMultiSet    = flags & (uint32)0x0001;        //  This is new in v02.

    _numFiles      = (uint64)1 << _numFilesBits;    //  The same for all formats, but
    _numBlocks     = (uint64)1 << _numBlocksBits;   //  awkward to do outside of here.
  }

  //  If we didn't initialize, set the file position to the start
  //  of the statistics.
  else {
    masterIndex->setPosition(64 + 64 + 32 + 32 + 32 + 32 + 32);
  }
}



void
kmerCountFileReader::initializeFromMasterI_v03(stuffedBits  *masterIndex,
                                               bool          doInitialize) {
  initializeFromMasterI_v02(masterIndex, doInitialize);
}



void
kmerCountFileReader::initializeFromMasterIndex(bool  doInitialize,
                                               bool  loadStatistics,
                                               bool  beVerbose) {
  char   N[FILENAME_MAX+1];

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _inName);

  if (fileExists(N) == false)
    fprintf(stderr, "ERROR: '%s' doesn't appear to be a meryl input; file '%s' doesn't exist.\n",
            _inName, N), exit(1);

  //  Open the master index.

  stuffedBits  *masterIndex = new stuffedBits(N);

  //  Based on the magic number, initialzie.

  uint64  m1 = masterIndex->getBinary(64);
  uint64  m2 = masterIndex->getBinary(64);
  uint32  vv = 1;

  if        ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
             (m2 == 0x31302e765f5f7865llu)) {   //  ex__v.01
    initializeFromMasterI_v01(masterIndex, doInitialize);
    vv = 1;

  } else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
             (m2 == 0x32302e765f5f7865llu)) {   //  ex__v.02
    initializeFromMasterI_v02(masterIndex, doInitialize);
    vv = 2;

  } else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
             (m2 == 0x33302e765f5f7865llu)) {   //  ex__v.03
    initializeFromMasterI_v03(masterIndex, doInitialize);
    vv = 3;

  } else {
    fprintf(stderr, "ERROR: '%s' doesn't look like a meryl input; file '%s' fails magic number check.\n",
            _inName, N), exit(1);
  }

  //  Check that the mersize is set and valid.

  uint32  merSize = (_prefixSize + _suffixSize) / 2;

  if (kmer::merSize() == 0)         //  If the global kmer size isn't set yet,
    kmer::setSize(merSize);         //  set it.

  if (kmer::merSize() != merSize)   //  And if set, make sure we're compatible.
    fprintf(stderr, "mer size mismatch, can't process this set of files.\n"), exit(1);

  //  If loading statistics is enabled, load the stats assuming the file is in
  //  the proper position.

  if (loadStatistics == true) {
    _stats = new kmerCountStatistics;
    _stats->load(masterIndex, vv);
  }

  //  And report some logging.

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

  delete masterIndex;
}



kmerCountFileReader::kmerCountFileReader(const char *inputName,
                                         bool        beVerbose) {
  strncpy(_inName, inputName, FILENAME_MAX);
  initializeFromMasterIndex(true, false, beVerbose);
}



kmerCountFileReader::kmerCountFileReader(const char *inputName,
                                         uint32      threadFile,
                                         bool        beVerbose) {
  strncpy(_inName, inputName, FILENAME_MAX);
  initializeFromMasterIndex(true, false, beVerbose);
  enableThreads(threadFile);
}



kmerCountFileReader::~kmerCountFileReader() {

  delete [] _blockIndex;

  delete [] _suffixes;
  delete [] _values;

  delete    _stats;

  AS_UTL_closeFile(_datFile);

  delete    _block;
}



void
kmerCountFileReader::loadStatistics(void) {
  if (_stats == NULL)
    initializeFromMasterIndex(false, true, false);
}



void
kmerCountFileReader::dropStatistics(void) {
  delete _stats;
  _stats = NULL;
}



void
kmerCountFileReader::enableThreads(uint32 threadFile) {
  _activeFile = threadFile;
  _threadFile = threadFile;
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
    _value = _values[_activeMer];
    return(true);
  }

  //  If no file, open whatever is 'active'.  In thread mode, the first file
  //  we open is the 'threadFile'; in normal mode, the first file we open is
  //  the first file in the database.

 loadAgain:
  if (_datFile == NULL)
    _datFile = openInputBlock(_inName, _activeFile, _numFiles);

  //  Load blocks.

  bool loaded = _block->loadBlock(_datFile, _activeFile);

  //  If nothing loaded. open a new file and try again.

  if (loaded == false) {
    AS_UTL_closeFile(_datFile);

    if (_activeFile == _threadFile)   //  Thread mode, if no block was loaded,
      return(false);                  //  we're done.

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

  resizeArrayPair(_suffixes, _values, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

  //  Decode the block into _OUR_ space.
  //
  //  decodeBlock() marks the block as having no data, so the next time we loadBlock() it will
  //  read more data from disk.  For blocks that don't get decoded, they retain whatever was
  //  loaded, and do not load another block in loadBlock().

  _block->decodeBlock(_suffixes, _values);

  //  But if no kmers in this block, load another block.  Sadly, the block must always
  //  be decoded, otherwise, the load will not load a new block.

  if (_nKmers == 0)
    goto loadAgain;

  //  Reset iteration, and load the first kmer.

  _activeMer = 0;

  _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
  _value = _values[_activeMer];

  return(true);
}
