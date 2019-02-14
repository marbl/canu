
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

//  Responsible for taking a block of sorted kmers
//  and writing them to the appropriate file.

class kmerCountFileWriter;

class kmerCountBlockWriter {
public:
  kmerCountBlockWriter(kmerCountFileWriter *writer);
  ~kmerCountBlockWriter();

public:
  void    addBlock(uint64  prefix, uint64  nKmers, uint64 *suffixes, uint32 *values);
  void    addBlock(uint64  prefix, uint64  nKmers, uint64 *suffixes, uint64 *values);

  void    finishBatch(void);
  void    finish(void);

private:
  void    closeFileDumpIndex(uint32 oi, uint32 iteration=UINT32_MAX);
  void    mergeBatches(uint32 oi);

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

  FILE                     **_datFiles;
  kmerCountFileIndex       **_datFileIndex;

  //  Kmer data and etc

  uint32                     _iteration;
};



