
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



void
kmerCountFileWriter::initialize(uint32 prefixSize, bool isMultiSet) {

  if (_initialized == true)    //  Nothing to do if we're already done.
    return;

  //  If the global mersize isn't set, we're hosed.

  if (kmer::merSize() == 0)
    fprintf(stderr, "kmerCountFileWriter::initialize()-- asked to initialize, but kmer::merSize() is zero!\n"), exit(1);

  //  The count operations call initialize() exactly once, but nextMer() calls
  //  it once per file and so we need some kind of concurrency control here.

#pragma omp critical (kmerCountFileWriterInit)
  if (_initialized == false) {

    //  If the prefixSize is zero, set it to (arbitrary) 1/4 the kmer size.
    //  This happens in the streaming writer (which is used when meryl does
    //  any non-count operation).  The prefixSize here just controls how
    //  often we dump blocks to the file.

    if (_prefixSize == 0)
      _prefixSize = prefixSize;

#warning how to set prefix size for streaming operations?
    if (_prefixSize == 0)
      _prefixSize = 12;  //max((uint32)8, 2 * kmer::merSize() / 3);

    _suffixSize         = 2 * kmer::merSize() - _prefixSize;
    _suffixMask         = uint64MASK(_suffixSize);

    //  Decide how many files to write.  We can make up to 2^32 files, but will
    //  run out of file handles _well_ before that.  For now, limit to 2^6 = 64 files.

    _numFilesBits       = 6;  //(_prefixSize < 7) ? _prefixSize : 6;
    _numBlocksBits      = _prefixSize - _numFilesBits;

    _numFiles           = (uint64)1 << _numFilesBits;
    _numBlocks          = (uint64)1 << _numBlocksBits;

    _isMultiSet         = isMultiSet;

    //  Now we're initialized!

    fprintf(stderr, "kmerCountFileWriter()-- Creating '%s' for %u-mers, with prefixSize %u suffixSize %u numFiles %lu\n",
            _outName, (_prefixSize + _suffixSize) / 2, _prefixSize, _suffixSize, _numFiles);

    _initialized = true;
  }
}



kmerCountFileWriter::kmerCountFileWriter(const char *outputName,
                                         uint32      prefixSize) {

  //  Note that we're not really initialized yet.  We could call initialize() in some cases,
  //  but the interesting one can't initialized() until the first meryl input file is opened,
  //  so we don't initialize any of them.

  _initialized   = false;

  //  Save the output directory name, and try to make it.  If we can't we'll fail quickly.

  strncpy(_outName, outputName, FILENAME_MAX);

  AS_UTL_mkdir(_outName);

  //  Parameters on how the suffixes/values are encoded are set once we know
  //  the kmer size.  See initialize().

  _prefixSize    = prefixSize;

  _suffixSize    = 0;
  _suffixMask    = 0;

  _numFilesBits  = 0;
  _numBlocksBits = 0;
  _numFiles      = 0;
  _numBlocks     = 0;

  _isMultiSet    = false;
}



kmerCountFileWriter::~kmerCountFileWriter() {
  uint32   flags = (uint32)0x0000;

  //  Set flags.

  if (_isMultiSet)
    flags |= (uint32)0x0001;

  //  Create a master index with the parameters.

  stuffedBits  *masterIndex = new stuffedBits;

  masterIndex->setBinary(64, 0x646e496c7972656dllu);  //  HEX: ........  ONDISK: merylInd
  masterIndex->setBinary(64, 0x33302e765f5f7865llu);  //       30.v__xe          ex__v.03
  masterIndex->setBinary(32, _prefixSize);
  masterIndex->setBinary(32, _suffixSize);
  masterIndex->setBinary(32, _numFilesBits);
  masterIndex->setBinary(32, _numBlocksBits);
  masterIndex->setBinary(32, flags);

  _stats.dump(masterIndex);

  //  Store the master index (and stats) to disk.

  char     N[FILENAME_MAX+1];
  FILE    *F;

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _outName);

  F = AS_UTL_openOutputFile(N);
  masterIndex->dumpToFile(F);
  AS_UTL_closeFile(F);

  delete masterIndex;
}



uint32
kmerCountFileWriter::fileNumber(uint64  prefix) {

  assert(_initialized);

  //  Based on the prefix, decide what output file to write to.
  //  The prefix has _prefixSize bits.  We want to save the highest _numFiles bits.

  uint64  oi  = prefix >> _numBlocksBits;

  if (oi >= _numFiles) {
    fprintf(stderr, "kmerCountFileWriter()-- Formed invalid file number %lu >= number of files %lu:\n", oi, _numFiles);
    fprintf(stderr, "kmerCountFileWriter()--   prefix          0x%016lx\n", prefix);
    fprintf(stderr, "kmerCountFileWriter()--   prefixSize      %u\n", _prefixSize);
    fprintf(stderr, "kmerCountFileWriter()--   suffixSize      %u\n", _suffixSize);
    fprintf(stderr, "kmerCountFileWriter()--   numFilesBits    %u\n", _numFilesBits);
    fprintf(stderr, "kmerCountFileWriter()--   numBlocksBits   %u\n", _numBlocksBits);
  }
  assert(oi < _numFiles);

  return((uint32)oi);
}



//void
//kmerCountFileWriter::importStatistics(kmerCountStatistics &import) {
//#warning NOT IMPORTING STATISTICS
//}



void
kmerCountFileWriter::writeBlockToFile(FILE                *datFile,
                                      kmerCountFileIndex  *datFileIndex,
                                      uint64               prefix,
                                      uint64               nKmers,
                                      uint64              *suffixes,
                                      uint32              *values) {

  //  Figure out the optimal size of the Elias-Fano prefix.  It's just log2(N)-1.

  uint32  unaryBits = 0;
  uint64  unarySum  = 1;
  while (unarySum < nKmers) {
    unaryBits  += 1;
    unarySum  <<= 1;
  }

  uint32  binaryBits = _suffixSize - unaryBits;      //  Only _suffixSize is used from the class.

  //  Dump data.

  stuffedBits   *dumpData = new stuffedBits;

  dumpData->setBinary(64, 0x7461446c7972656dllu);    //  Magic number, part 1.
  dumpData->setBinary(64, 0x0a3030656c694661llu);    //  Magic number, part 2.

  dumpData->setBinary(64, prefix);
  dumpData->setBinary(64, nKmers);

  dumpData->setBinary(8,  1);                        //  Kmer coding type
  dumpData->setBinary(32, unaryBits);                //  Kmer coding parameters
  dumpData->setBinary(32, binaryBits);
  dumpData->setBinary(64, 0);

  dumpData->setBinary(8,  1);                        //  Value coding type
  dumpData->setBinary(64, 0);                        //  Value coding parameters
  dumpData->setBinary(64, 0);

  //  Split the kmer suffix into two pieces, one unary encoded offsets and one binary encoded.

  uint64  lastPrefix = 0;
  uint64  thisPrefix = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    thisPrefix = suffixes[kk] >> binaryBits;

    dumpData->setUnary(thisPrefix - lastPrefix);
    dumpData->setBinary(binaryBits, suffixes[kk]);

    lastPrefix = thisPrefix;
  }

  //  Save the values, too.  Eventually these will be cleverly encoded.  Really.

  uint64  lastValue = 0;
  uint64  thisValue = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    dumpData->setBinary(32, values[kk]);
  }

  //  Save the index entry.

  uint64  block = prefix & uint64MASK(_numBlocksBits);

  datFileIndex[block].set(prefix, datFile, nKmers);

  //  Dump data to disk, cleanup, and done!

  dumpData->dumpToFile(datFile);

  delete dumpData;
}



void
kmerCountFileWriter::writeBlockToFile(FILE                *datFile,
                                      kmerCountFileIndex  *datFileIndex,
                                      uint64               prefix,
                                      uint64               nKmers,
                                      uint64              *suffixes,
                                      uint64              *values) {

  //  Figure out the optimal size of the Elias-Fano prefix.  It's just log2(N)-1.

  uint32  unaryBits = 0;
  uint64  unarySum  = 1;
  while (unarySum < nKmers) {
    unaryBits  += 1;
    unarySum  <<= 1;
  }

  uint32  binaryBits = _suffixSize - unaryBits;      //  Only _suffixSize is used from the class.

  //  Dump data.

  stuffedBits   *dumpData = new stuffedBits;

  dumpData->setBinary(64, 0x7461446c7972656dllu);    //  Magic number, part 1.
  dumpData->setBinary(64, 0x0a3030656c694661llu);    //  Magic number, part 2.

  dumpData->setBinary(64, prefix);
  dumpData->setBinary(64, nKmers);

  dumpData->setBinary(8,  1);                        //  Kmer coding type
  dumpData->setBinary(32, unaryBits);                //  Kmer coding parameters
  dumpData->setBinary(32, binaryBits);
  dumpData->setBinary(64, 0);

  dumpData->setBinary(8,  2);                        //  Value coding type
  dumpData->setBinary(64, 0);                        //  Value coding parameters
  dumpData->setBinary(64, 0);

  //  Split the kmer suffix into two pieces, one unary encoded offsets and one binary encoded.

  uint64  lastPrefix = 0;
  uint64  thisPrefix = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    thisPrefix = suffixes[kk] >> binaryBits;

    dumpData->setUnary(thisPrefix - lastPrefix);
    dumpData->setBinary(binaryBits, suffixes[kk]);

    lastPrefix = thisPrefix;
  }

  //  Save the values, too.  Eventually these will be cleverly encoded.  Really.

  uint64  lastValue = 0;
  uint64  thisValue = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    dumpData->setBinary(64, values[kk]);
  }

  //  Save the index entry.

  uint64  block = prefix & uint64MASK(_numBlocksBits);

  datFileIndex[block].set(prefix, datFile, nKmers);

  //  Dump data to disk, cleanup, and done!

  dumpData->dumpToFile(datFile);

  delete dumpData;
}
