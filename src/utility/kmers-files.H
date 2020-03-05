
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

#ifndef LIBKMER_FILES_H
#define LIBKMER_FILES_H

#ifndef LIBKMER_H
#error "include kmers.H, not this."
#endif



void
dumpMerylDataFile(char *name);



//  Functions to constrct data file names and open them for reading or
//  writing.

char *
constructBlockName(char   *nameprefix,
                   uint64  outIndex,
                   uint32  numFiles,
                   uint32  iteration,
                   bool    isIndex);

FILE *
openOutputBlock(char   *nameprefix,
                uint64  fileIndex,
                uint32  numFiles,
                uint32  iteration=0);

FILE *
openInputBlock(char   *nameprefix,
               uint64  fileIndex,
               uint32  numFiles,
               uint32  iteration=0);


//  Read a block of kmer data from disk, and decode it into a list of kmers,
//  counts and (eventually) colors.

class merylFileBlockReader {
public:
  merylFileBlockReader();
  ~merylFileBlockReader();

  bool      loadBlock(FILE *inFile, uint32 activeFile, uint32 activeIteration=0);

  void      decodeBlock(void);                               //  to our own storage
  void      decodeBlock(kmdata *suffixes, kmvalu *values);   //  to external storage

  kmpref    prefix(void)   { return(_blockPrefix); };        //  kmer prefix of this block
  uint64    nKmers(void)   { return(_nKmers);      };        //  number of kmers in this block

  kmdata   *suffixes(void) { return(_suffixes); };           //  direct access to decoded data
  kmvalu   *values(void)   { return(_values);   };

private:
  stuffedBits  *_data;

  kmpref        _blockPrefix;  //  The prefix of all kmers in this block
  uint64        _nKmers;       //  The number of kmers in this block
  uint64        _nKmersMax;    //  The number of kmers we've allocated space for in _suffixes and _values

  uint32        _kCode;        //  Encoding type of kmer, then 128 bits of parameters
  uint32        _unaryBits;    //    bits in the unary prefix  (of the kmer suffix)
  uint32        _binaryBits;   //    bits in the binary suffix (of the kmer suffix)
  uint64        _k1;           //    unused

  uint32        _cCode;        //  Encoding type of the values, then 128 bits of parameters
  uint64        _c1;           //    unused
  uint64        _c2;           //    unused

  kmdata       *_suffixes;     //  Decoded suffixes and values.
  kmvalu       *_values;       //
};


//  An index to the binary encoded kmer data.  Provides:
//    kmer prefix for each block
//    starting position of the block in the file
//    number of kmers in each block
//
//  Used as argument to merylFileReader::loadFromFile().
//  Populated by the file writer.

class merylFileIndex {
public:
  merylFileIndex() {
    clear();
  };

  void       set(kmpref  prefix,
                 FILE   *F,
                 uint64  nKmers) {

    if (_blockPosition == UINT64_MAX) {
      _blockPrefix   = prefix;
      _blockPosition = AS_UTL_ftell(F);
      _numKmers      = nKmers;
    }

    else {
      _numKmers     += nKmers;
    }

    if (_blockPrefix != prefix)
      fprintf(stderr, "set prefix 0x%s from [0x%s %8lu %8lu] to [%8lu %8lu]\n",
              toHex(prefix),
              toHex(_blockPrefix), _blockPosition, _numKmers,
              AS_UTL_ftell(F), _numKmers + nKmers);
    assert(_blockPrefix   == prefix);

    assert(_blockPosition <= AS_UTL_ftell(F));
  };

  void       clear(void) {
    _blockPrefix   = 0;
    _blockPosition = UINT64_MAX;
    _numKmers      = 0;
  }

  kmpref     blockPrefix(void)     { return((kmpref)_blockPrefix);   };
  uint64     blockPosition(void)   { return(        _blockPosition); };
  uint64     numKmers(void)        { return(        _numKmers);      };

private:
  uint64    _blockPrefix;     //  For compatibility, and alignment, _blockPrefix
  uint64    _blockPosition;   //  needs to be uint64 instead of the more correct
  uint64    _numKmers;        //  kmpref.
};



#endif  //  LIBKMER_FILES_H
