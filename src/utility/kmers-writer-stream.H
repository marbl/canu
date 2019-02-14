
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

//  Responsible for taking one kmer at a time and writing
//  it to

class kmerCountFileWriter;

class kmerCountStreamWriter {
public:
  kmerCountStreamWriter(kmerCountFileWriter *writer,
                        uint32               fileNumber);

  ~kmerCountStreamWriter();

public:
  void    addMer(kmer k, uint64 c);

private:
  void    dumpBlock(uint64 nextPrefix=UINT64_MAX);

private:
  kmerCountFileWriter       *_writer;
  char                       _outName[FILENAME_MAX+1];

  //  Encoding data

  uint32                     _prefixSize;

  uint32                     _suffixSize;
  uint64                     _suffixMask;

  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;
  uint64                     _numFiles;
  uint64                     _numBlocks;

  //  File data

  uint32                     _filePrefix;

  FILE                      *_datFile;
  kmerCountFileIndex        *_datFileIndex;

  //  Kmer data and etc for writing the stream.

  uint64                     _batchPrefix;
  uint64                     _batchNumKmers;
  uint64                     _batchMaxKmers;
  uint64                    *_batchSuffixes;
  uint64                    *_batchValues;

  //kmerCountStatistics        _stats;
};



