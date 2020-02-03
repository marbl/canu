
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

#ifndef LIBKMER_FILES_READER_H
#define LIBKMER_FILES_READER_H

#ifndef LIBKMER_H
#error "include kmers.H, not this."
#endif


class merylFileReader {
private:
  void    initializeFromMasterI_v00(void);
  void    initializeFromMasterI_v01(stuffedBits  *masterIndex, bool doInitialize);
  void    initializeFromMasterI_v02(stuffedBits  *masterIndex, bool doInitialize);
  void    initializeFromMasterI_v03(stuffedBits  *masterIndex, bool doInitialize);
  void    initializeFromMasterIndex(bool  doInitialize, bool  loadStatistics, bool  beVerbose);

public:
  merylFileReader(const char *inputName,
                      bool        beVerbose=false);

  merylFileReader(const char *inputName,
                      uint32      threadFile,
                      bool        beVerbose=false);

  ~merylFileReader();

public:
  void    rewind(void) {
    _activeMer  = 0;    //  Position we are at in the block loaded.
    _activeFile = 0;

    _nKmers     = 0;    //  Kmers in the block we have loaded.

    if (_threadFile != UINT32_MAX)
      _activeFile = _threadFile;

    AS_UTL_closeFile(_datFile);
  };

public:
  void    loadStatistics(void);
  void    dropStatistics(void);

public:
  void    enableThreads(uint32 threadFile);

public:
  void    loadBlockIndex(void);

public:
  bool    nextMer(void);
  kmer    theFMer(void)        { return(_kmer);        };
  kmvalu  theValue(void)       { return(_value);       };

  bool    isMultiSet(void)     { return(_isMultiSet);  };

  char   *filename(void)       { return(_inName);      };

  merylHistogram  *stats(void) {
    loadStatistics();
    return(_stats);
  }

  //  For direct access to the kmer blocks.
public:
  uint32  prefixSize(void)     { return(_prefixSize); };
  uint32  suffixSize(void)     { return(_suffixSize); };

  uint32  numFiles(void)       { return(_numFiles);  };
  uint32  numBlocks(void)      { return(_numBlocks); };

  uint32  numFilesBits(void)   { return(_numFilesBits);  };
  uint32  numBlocksBits(void)  { return(_numBlocksBits); };

  FILE   *blockFile(uint32 ff) {
    FILE *F = NULL;

    if (ff < _numFiles)
      F = openInputBlock(_inName, ff, _numFiles);

    return(F);
  };

  merylFileIndex   &blockIndex(uint32 bb) {
    return(_blockIndex[bb]);
  };

private:
  char                       _inName[FILENAME_MAX+1];

  uint32                     _prefixSize;
  uint32                     _suffixSize;
  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;

  uint32                     _numFiles;
  uint32                     _numBlocks;

  bool                       _isMultiSet;

  merylHistogram            *_stats;

  FILE                      *_datFile;

  merylFileBlockReader      *_block;
  merylFileIndex            *_blockIndex;

  kmer                       _kmer;
  kmvalu                     _value;

  uint64                     _prefix;

  kmdata                     _activeMer;
  uint32                     _activeFile;

  uint32                     _threadFile;

  uint64                     _nKmers;
  uint64                     _nKmersMax;
  kmdata                    *_suffixes;
  kmvalu                    *_values;
};



#endif  //  LIBKMER_FILES_READER_H
