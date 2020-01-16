
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
 *    Brian P. Walenz beginning on 2018-AUG-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kmers.H"
#include "bits.H"

#include "files.H"



merylStreamWriter::merylStreamWriter(merylFileWriter *writer, uint32 fileNumber) {

  _writer = writer;

  strncpy(_outName, _writer->_outName, FILENAME_MAX);

  //  Encoding data

  _prefixSize    = _writer->_prefixSize;

  _suffixSize    = _writer->_suffixSize;
  _suffixMask    = _writer->_suffixMask;

  _numFilesBits  = _writer->_numFilesBits;
  _numBlocksBits = _writer->_numBlocksBits;
  _numFiles      = _writer->_numFiles;
  _numBlocks     = _writer->_numBlocks;

  //  File data

  _filePrefix    = fileNumber;

  _datFile       = openOutputBlock(_outName, _filePrefix, _numFiles, 0);
  _datFileIndex  = new merylFileIndex [_numBlocks];

  //  Kmer data

  _batchPrefix   = 0;
  _batchNumKmers = 0;
  _batchMaxKmers = 16 * 1048576;
  _batchSuffixes = NULL;
  _batchValues   = NULL;
}



merylStreamWriter::~merylStreamWriter() {

  //  If data in the batch, dump it.  Cleanup and close the data file.

  if (_batchNumKmers > 0)
    dumpBlock();

  delete [] _batchSuffixes;
  delete [] _batchValues;

  AS_UTL_closeFile(_datFile);

  //  Write the index data for this file.

  char  *idxname = constructBlockName(_outName, _filePrefix, _numFiles, 0, true);
  FILE  *idxfile = AS_UTL_openOutputFile(idxname);

  writeToFile(_datFileIndex, "merylStreamWriter::fileIndex", _numBlocks, idxfile);

  AS_UTL_closeFile(idxfile, idxname);

  delete [] idxname;

  delete [] _datFileIndex;
}



void
merylStreamWriter::dumpBlock(kmpref nextPrefix) {

  //fprintf(stderr, "merylStreamWriter::dumpBlock()-- write batch for prefix %lu with %lu kmers.\n",
  //        _batchPrefix, _batchNumKmers);

  //  Encode and dump to disk.

  _writer->writeBlockToFile(_datFile, _datFileIndex,
                            _batchPrefix,
                            _batchNumKmers,
                            _batchSuffixes,
                            _batchValues);

  //  Insert counts into the histogram.

#pragma omp critical (merylFileWriterAddValue)
  for (uint32 kk=0; kk<_batchNumKmers; kk++)
    _writer->_stats.addValue(_batchValues[kk]);

  //  Set up for the next block of kmers.

  _batchPrefix   = nextPrefix;
  _batchNumKmers = 0;
}



void
merylStreamWriter::addMer(kmer k, kmvalu c) {

  kmpref  prefix = (kmdata)k >> _suffixSize;
  kmdata  suffix = (kmdata)k  & _suffixMask;

  //  Do we need to initialize to firstPrefixInFile(ff) and also write empty prefixes?
  //  Or can we just init to the first prefix we see?

  if (_batchSuffixes == NULL) {
    //fprintf(stderr, "merylFileWriter::addMer()-- ff %2u allocate %7lu kmers for a batch\n", ff, _batchMaxKmers);
    _batchPrefix   = prefix;
    _batchNumKmers = 0;
    _batchMaxKmers = 16 * 1048576;
    _batchSuffixes = new kmdata [_batchMaxKmers];
    _batchValues   = new kmvalu [_batchMaxKmers];
  }

  //  If the batch is full, or we've got a kmer for a different batch, dump the batch
  //  to disk.

  bool  dump1 = (_batchNumKmers >= _batchMaxKmers);
  bool  dump2 = (_batchPrefix != prefix) && (_batchNumKmers > 0);

  //if (dump1 || dump2)
  //  fprintf(stderr, "merylStreamWriter::addMer()-- ff %u addBlock 0x%016lx with %lu kmers dump1=%c dump2=%c\n",
  //          ff, _batchPrefix, _batchNumKmers,
  //          (dump1) ? 'T' : 'F',
  //          (dump2) ? 'T' : 'F');

  if (dump1 || dump2)
    dumpBlock(prefix);

  //  And now just add the kmer to the list.

  assert(_batchNumKmers < _batchMaxKmers);

  _batchSuffixes[_batchNumKmers] = suffix;
  _batchValues  [_batchNumKmers] = c;

  _batchNumKmers++;
}
