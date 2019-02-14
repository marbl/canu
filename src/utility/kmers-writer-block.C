
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



kmerCountBlockWriter::kmerCountBlockWriter(kmerCountFileWriter *writer) {

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

  _datFiles      = new FILE *               [_numFiles];
  _datFileIndex  = new kmerCountFileIndex * [_numFiles];

  for (uint32 ii=0; ii<_numFiles; ii++) {
    _datFiles[ii]     = NULL;
    _datFileIndex[ii] = new kmerCountFileIndex [_numBlocks];
  }

  //  Kmer data

  _iteration     = 1;
}



kmerCountBlockWriter::~kmerCountBlockWriter() {

  //  Check that all output files are closed.  If not,
  //  finishIteration() was never called.  Since that's somewhat heavy-weight,
  //  and threaded, we don't want to call it now.

  for (uint32 ii=0; ii<_numFiles; ii++) {
    if (_datFiles[ii] != NULL)
      fprintf(stderr, "ERROR: open output file %u exists when destroying writer for '%s'\n",
              ii, _outName);
    assert(_datFiles[ii] == NULL);
  }

  delete [] _datFiles;

  //  Write index data.
#if 0
  for (uint32 ii=0; ii<_numFiles; ii++) {
    char  *idxname = constructBlockName(_outName, ii, _numFiles, 0, true);
    FILE  *idxfile = AS_UTL_openOutputFile(idxname);

    writeToFile(_datFileIndex[ii], "kmerCountBlockWriter::fileIndex", _numBlocks, idxfile);

    AS_UTL_closeFile(idxfile, idxname);

    delete [] idxname;
  }
#endif

  for (uint32 ii=0; ii<_numFiles; ii++)
    delete [] _datFileIndex[ii];

  delete [] _datFileIndex;

  //  Tell the master that we're done.

  //for (uint32 ii=0; ii<_numFiles; ii++)
  //  _writer->_stats.import(_stats);
}



//  We're given a block of kmers with the same prefix, sorted, with values.
//  Which is exactly what kmerCountFileWriter wants, so we can just forward
//  the call to it.
//
void
kmerCountBlockWriter::addBlock(uint64  prefix,
                               uint64  nKmers,
                               uint64 *suffixes,
                               uint32 *values) {

  //  Open a new file, if needed.

  uint32 oi = _writer->fileNumber(prefix);

  if (_datFiles[oi] == NULL)
    _datFiles[oi] = openOutputBlock(_outName, oi, _numFiles, _iteration);

  //  Encode and dump to disk.

  _writer->writeBlockToFile(_datFiles[oi], _datFileIndex[oi],
                            prefix,
                            nKmers,
                            suffixes,
                            values);

  //  Insert values into the histogram.

#pragma omp critical (kmerCountFileWriterAddValue)
  for (uint32 kk=0; kk<nKmers; kk++)
    _writer->_stats.addValue(values[kk]);
}



void
kmerCountBlockWriter::addBlock(uint64  prefix,
                               uint64  nKmers,
                               uint64 *suffixes,
                               uint64 *values) {

  //  Open a new file, if needed.

  uint32 oi = _writer->fileNumber(prefix);

  if (_datFiles[oi] == NULL)
    _datFiles[oi] = openOutputBlock(_outName, oi, _numFiles, _iteration);

  //  Encode and dump to disk.

  _writer->writeBlockToFile(_datFiles[oi], _datFileIndex[oi],
                            prefix,
                            nKmers,
                            suffixes,
                            values);

  //  Insert values into the histogram.

#pragma omp critical (kmerCountFileWriterAddValue)
  for (uint32 kk=0; kk<nKmers; kk++)
    _writer->_stats.addValue(values[kk]);
}



//  Close all data files.
//  Write their merylIndex files.
//  Clear the index data in memory.
void
kmerCountBlockWriter::closeFileDumpIndex(uint32 oi, uint32 iteration) {

  //  Close all data files.

  AS_UTL_closeFile(_datFiles[oi]);

  //  Write block indices, then clear each one.

  if (iteration == UINT32_MAX)
    iteration = _iteration;

  char  *idxname = constructBlockName(_outName, oi, _numFiles, iteration, true);
  FILE  *idxfile = AS_UTL_openOutputFile(idxname);

  writeToFile(_datFileIndex[oi], "kmerCountBlockWriter::fileIndex", _numBlocks, idxfile);

  AS_UTL_closeFile(idxfile, idxname);

  delete [] idxname;

  for (uint32 bb=0; bb<_numBlocks; bb++)
    _datFileIndex[oi][bb].clear();
}



void
kmerCountBlockWriter::finishBatch(void) {

  for (uint32 ii=0; ii<_numFiles; ii++)
    closeFileDumpIndex(ii);

  _iteration++;
}



void
kmerCountBlockWriter::finish(void) {

  fprintf(stderr, "finishIteration()--\n");

  for (uint32 ii=0; ii<_numFiles; ii++)
    closeFileDumpIndex(ii);

  //  If only one iteration, just rename files to the proper name.

  if (_iteration == 1) {
    char *oldName;
    char *newName;

    for (uint32 oi=0; oi<_numFiles; oi++) {
      oldName = constructBlockName(_outName, oi, _numFiles, 1, false);  //  Data files.
      newName = constructBlockName(_outName, oi, _numFiles, 0, false);

      AS_UTL_rename(oldName, newName);

      delete [] newName;
      delete [] oldName;

      oldName = constructBlockName(_outName, oi, _numFiles, 1, true);  //  Index files.
      newName = constructBlockName(_outName, oi, _numFiles, 0, true);

      AS_UTL_rename(oldName, newName);

      delete [] newName;
      delete [] oldName;
    }
  }

  //  Otherwise, merge the multiple iterations into a single file (after clearing
  //  the (unused) stats from the last iteration).

  else {
    _writer->_stats.clear();

    fprintf(stderr, "finishIteration()--  Merging %u blocks.\n", _iteration);

#pragma omp parallel for
    for (uint32 oi=0; oi<_numFiles; oi++)
      mergeBatches(oi);
  }
}



void
kmerCountBlockWriter::mergeBatches(uint32 oi) {
  kmerCountFileReaderBlock    inBlocks[_iteration + 1];
  FILE                       *inFiles [_iteration + 1];

  //  Open the input files, allocate blocks.

  inFiles[0] = NULL;

  for (uint32 ii=1; ii <= _iteration; ii++)
    inFiles[ii]  = openInputBlock(_outName, oi, _numFiles, ii);

  //  Open the output file.

  assert(_datFiles[oi] == NULL);

  _datFiles[oi] = openOutputBlock(_outName, oi, _numFiles);

  //  Create space to save out suffixes and values.

  uint64    nKmersMax = 0;
  uint64   *suffixes  = NULL;
  uint64   *values    = NULL;

  uint64    kmersIn   = 0;
  uint64    kmersOut  = 0;

  //  Load each block from each file, merge, and write.

  uint32    p[_iteration+1];  //  Position in s[] and v[]
  uint64    l[_iteration+1];  //  Number of entries in s[] and v[]
  uint64   *s[_iteration+1];  //  Pointer to the suffixes for piece x
  uint64   *v[_iteration+1];  //  Pointer to the values   for piece x

  for (uint32 bb=0; bb<_numBlocks; bb++) {
    uint64  totnKmers = 0;
    uint64  savnKmers = 0;

    //  Load and decode each block.  NO ERROR CHECKING.

    for (uint32 ii=1; ii <= _iteration; ii++) {        //  This loop _could_ be threaded, but the
      inBlocks[ii].loadBlock(inFiles[ii], oi, ii);     //  caller to this function has already
      inBlocks[ii].decodeBlock();                      //  threaded things, so no real point.

      p[ii] = 0;
      l[ii] = inBlocks[ii].nKmers();
      s[ii] = inBlocks[ii].suffixes();
      v[ii] = inBlocks[ii].values();

      totnKmers += l[ii];
    }

    //  Check that everyone has loaded the same prefix.

    uint64  prefix = inBlocks[1].prefix();

    //fprintf(stderr, "MERGE prefix 0x%04lx %8lu kmers and 0x%04lx %8lu kmers.\n",
    //        prefix, l[1], inBlocks[2].prefix(), l[2]);

    for (uint32 ii=1; ii <= _iteration; ii++) {
      if (prefix != inBlocks[ii].prefix())
        fprintf(stderr, "ERROR: File %u segments 1 and %u differ in prefix: 0x%016lx vs 0x%016lx\n",
                oi, ii, prefix, inBlocks[ii].prefix());
      assert(prefix == inBlocks[ii].prefix());
    }

    //  Setup the merge.

    resizeArrayPair(suffixes, values, 0, nKmersMax, totnKmers, resizeArray_doNothing);

    //  Merge!  We don't know the number of different kmers in the input, and are forced
    //  to loop infinitely.

    while (1) {
      uint64  minSuffix = UINT64_MAX;
      uint64  sumValue  = 0;

      //  Find the smallest suffix over all the inputs;
      //  Remember the sum of their values.

      for (uint32 ii=1; ii <= _iteration; ii++) {
        if (p[ii] < l[ii]) {
          if (minSuffix > s[ii][ p[ii] ]) {
            minSuffix = s[ii][ p[ii] ];
            sumValue  = v[ii][ p[ii] ];
          }

          else if (minSuffix == s[ii][ p[ii] ]) {
            sumValue += v[ii][ p[ii] ];
          }
        }
      }

      //  If no values, we're done.

      if ((minSuffix == UINT64_MAX) && (sumValue == 0))
        break;

      //  Set the suffix/value in our merged list, reallocating if needed.

      suffixes[savnKmers] = minSuffix;
      values  [savnKmers] = sumValue;

      savnKmers++;

      if (savnKmers > nKmersMax)
        fprintf(stderr, "savnKmers %lu > nKmersMax %lu\n", savnKmers, nKmersMax);
      assert(savnKmers <= nKmersMax);

      //  Move to the next element of the lists we pulled data from.  If the list is
      //  exhausted, mark it as so.

      for (uint32 ii=1; ii <= _iteration; ii++) {
        if ((p[ii] < l[ii]) &&
            (minSuffix == s[ii][ p[ii] ]))
          p[ii]++;
      }
    }

    //  Write the merged block of data to the output.

    _writer->writeBlockToFile(_datFiles[oi], _datFileIndex[oi],
                              prefix,
                              savnKmers,
                              suffixes,
                              values);

    //  Finally, don't forget to insert the values into the histogram!

#pragma omp critical (kmerCountBlockWriterAddValue)
    for (uint32 kk=0; kk<savnKmers; kk++)
      _writer->_stats.addValue(values[kk]);

    //  And update our local stats

    kmersIn  += totnKmers;
    kmersOut += savnKmers;
  }

  delete [] suffixes;
  delete [] values;

  //  Close the input data files.

  for (uint32 ii=1; ii <= _iteration; ii++)
    AS_UTL_closeFile(inFiles[ii]);

  //  Close the output files.

  closeFileDumpIndex(oi, 0);

  //  Remove the old data files.

  for (uint32 ii=1; ii <= _iteration; ii++) {
    char    *dname = constructBlockName(_outName, oi, _numFiles, ii, false);  //  Data files.
    char    *iname = constructBlockName(_outName, oi, _numFiles, ii, true);   //  Index files.

    AS_UTL_unlink(dname);
    AS_UTL_unlink(iname);

    delete [] dname;
    delete [] iname;
  }
}

