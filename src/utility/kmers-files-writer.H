
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

#ifndef LIBKMER_FILES_WRITER_H
#define LIBKMER_FILES_WRITER_H

#ifndef LIBKMER_H
#error "include kmers.H, not this."
#endif

#include "kmers-files.H"
#include "kmers-writer-block.H"
#include "kmers-writer-stream.H"

#include "kmers-histogram.H"

//  Responsible for writing data to a file.
//
//  It's a little confusing to use; the merylFileWriter is used to get a
//  blockWriter or streamWriter object that actually accepts the data to
//  write.

class merylFileWriter {
public:
  merylFileWriter(const char *outputName,
                      uint32      prefixSize = 0);

  ~merylFileWriter();

public:
  char   *filename(void)  { return(_outName);  };

public:
  void    initialize(uint32 prefixSize = 0, bool isMultiSet = false);

  merylBlockWriter  *getBlockWriter(void)        { return(new merylBlockWriter (this));      };
  merylStreamWriter *getStreamWriter(uint32 ff)  { return(new merylStreamWriter(this, ff));  };

public:
  uint32  numberOfFiles(void)           { return(_numFiles);                      };
  uint64  firstPrefixInFile(uint32 ff)  { return(((uint64)ff) << _numBlocksBits); };
  uint64  lastPrefixInFile(uint32 ff)   { return(firstPrefixInFile(ff + 1) - 1);  };

  uint32  fileNumber(uint64 prefix);

private:
  void    writeBlockToFile(FILE            *datFile,
                           merylFileIndex  *datFileIndex,
                           kmpref           blockPrefix,
                           uint64           nKmers,
                           kmdata          *suffixes,
                           uint32          *values);

  void    writeBlockToFile(FILE            *datFile,
                           merylFileIndex  *datFileIndex,
                           kmpref           blockPrefix,
                           uint64           nKmers,
                           kmdata          *suffixes,
                           uint64          *values);

private:
  bool                       _initialized;

  char                       _outName[FILENAME_MAX+1];

  uint32                     _prefixSize;

  uint32                     _suffixSize;
  uint64                     _suffixMask;

  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;
  uint64                     _numFiles;
  uint64                     _numBlocks;

  bool                       _isMultiSet;

  merylHistogram             _stats;

  friend class merylBlockWriter;
  friend class merylStreamWriter;
};




#endif  //  LIBKMER_FILES_WRITER_H
